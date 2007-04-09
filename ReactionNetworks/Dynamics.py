"""
Methods for evaluating the dynamics of Network.
"""
__docformat__ = "restructuredtext en"

import copy
import sets

import scipy
import scipy.optimize

import logging
logger = logging.getLogger('ReactionNetworks.Dynamics')

import SloppyCell.daskr as daskr
daeint = daskr.daeint
import Trajectory_mod

import SloppyCell.Utility as Utility
from SloppyCell import HAVE_PYPAR, my_rank, my_host, num_procs
if HAVE_PYPAR:
    import pypar

# XXX: Need to incorporate root_grace_t into integrate_sensitivity
global_rtol = 1e-6

return_derivs = False # do we want time derivatives of all trajectories returned?
reduce_space = 0 # we may want to not return every timepoint but only 0, reduce_space, 2*reduce_space etc.

failed_args = []
failed_kwargs = {}

class IntegrationException(Utility.SloppyCellException):
    pass

class FixedPointException(Utility.SloppyCellException):
    pass

def integrate_tidbit(net, res_func, Dfun, root_func, IC, yp0, curTimes, 
                     rtol, atol, fill_traj, return_derivs,
                     redirect_msgs, init_consistent, var_types):

    if root_func is not None:
        root_term = [False] * len(net._root_func)
        for ii in range(len(net.events)):
            root_term[ii] = True
    else:
        root_term = []

    N_dyn_var = len(net.dynamicVars)
    if net.integrateWithLogs:
        yp0 = yp0/IC
        IC = scipy.log(IC)

        old_func = res_func
        def res_func(time, logy, logyprime, ires):
            # Convert from logs of dyamic vars. We're integrating in logs
            #  partially to avoid zeros, so we put a lower bound on our
            #  output dynamicVars
            y = scipy.exp(logy)
            y = scipy.maximum(y, scipy.misc.limits.double_tiny)
            yprime = logyprime * y

            return old_func(time, y, yprime, ires)

        Dfun = None

        old_root_func = root_func
        def root_func(time, logy, logyprime):
            y = scipy.exp(logy)
            y = scipy.maximum(y, scipy.misc.limits.double_tiny)
            yprime = logyprime * y
            return old_root_func(time, y, yprime)

        atol = rtol
        rtol = [0]*len(rtol)

    int_args = {'res': res_func,
                't': curTimes,
                'y0': IC,
                'yp0': yp0,

                'jac': Dfun,
                'rt': root_func,

                'rtol': rtol,
                'atol': atol,

                'intermediate_output': fill_traj,
                'redir_output': redirect_msgs,

                'init_consistent' : init_consistent,
                'var_types' : var_types,
                }

    exception_raised = False
    try:
        int_returns = daeint(**int_args)
    except Utility.SloppyCellException, X:
        # When an exception happens, we'll try to return all the
        # trajectory we can, thus we don't reraise yet.
        exception_raised = X
        # We store the arguments to daeint that failed as
        #  Dynamics.failed_kwargs for debugging purproses.
        global failed_kwargs
        failed_kwargs = int_args
        # Recover what we have of int_returns from the exception,
        #  so we can process it.
        int_returns = X.args[1]

    # Pull things out of what daeint returned
    y_this, t_this = int_returns[0], int_returns[1]
    if return_derivs:
        ydt_this = int_returns[2]
    else:
        ydt_this = []

    t_root_this, y_root_this, i_root_this = int_returns[3:6]

    info_dict = int_returns[-1]

    # Deal with integrateWithLogs
    if net.integrateWithLogs:
        y_this = scipy.exp(y_this)
        y_root_this = scipy.exp(y_root_this)

    return exception_raised, y_this, t_this, ydt_this,\
            t_root_this, y_root_this, i_root_this

def generate_tolerances(net, rtol, atol=None):
    if rtol == None:
        rtol = 1e-6
    if scipy.isscalar(rtol):
        rtol = [rtol] * len(net.dynamicVars)
    rtol = scipy.minimum(rtol, global_rtol)

    # We set atol to be a minimum of 1e-6 to avoid variables with large typical
    #  values going negative.
    if (scipy.isscalar(atol) or atol==None):
        typ_vals = [net.get_var_typical_val(id) for id in net.dynamicVars.keys()]
        atol = scipy.minimum(rtol * scipy.asarray(typ_vals), 1e-6)
    return rtol, atol


def integrate(net, times, rtol=None, atol=None, params=None, fill_traj=True,
              return_events=False, return_derivs=False,
              redirect_msgs=True, calculate_ic = False,
              root_grace=True):
    """
    Integrate a Network, returning a Trajectory.

    net            The Network to integrate.
    times          A sequence of times to include in the integration output.
    params         Parameters for the integration. If params=None, the current
                   parameters of net are used.
    rtol           Relative error tolerance for this integration.  If rtol = None
                   then relative toleran
    atol           Absolute error tolerance for this integration.
    fill_traj      If True, the integrator will add the times it samples to
                   the returned trajectory. This is slower, but the resulting
                   trajectory will be more densely sampled where the dynamics
                   change quickly.
    return_events  If True, the output is (traj, te, ye, ie). te is a sequence
                   of times at which events fired. ye is a sequence of the
                   dynamic variable values where the events fired. ie is the
                   index of the fired events.
    return_derivs  If True, the returned trajectory records the time derivatives
                   of all quantities.
    redirect_msgs  If False, the errors and other output generated by the
                   integrator will be returned to the display.
    calculate_ic   If True, the integrator will calculate consistent initial
                   conditions.
    root_grace     If True, the integration won't fire events for a small period
                   after an event fires to prevent events from repeatedly firing
                   inapropriately.                   
                   
    """

    if params is not None:
        net.update_optimizable_vars(params)
    # If you ask for time = 0, we'll assume you want dynamic variable values
    # reset.
    if times[0] == 0:
        net.resetDynamicVariables()
    net.compile()

    if (rtol == None or atol == None):
        rtol, atol = generate_tolerances(net, rtol, atol)

    # if times isn't already an array, make it an array
    times = scipy.asarray(times)

    # get the initial state, time, and derivative.
    IC = net.getDynamicVarValues()
    start = times[0]
    # time derivative of trajectory
    youtdt = scipy.zeros((0, len(IC)), scipy.float_)
    ypIC = net.get_ddv_dt(IC,start)
    
    # start variables for output storage 
    yout = scipy.zeros((0, len(IC)), scipy.float_)
    tout = []

    # For some reason, the F2Py wrapper for ddaskr seems to fail if this
    #  isn't wrapped up in a lambda...
    res_func = lambda t, y, yprime, ires: net.res_function(t, y, yprime, ires)
    root_func = net.ddaskr_root
    num_events = len(net.events)
    
    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events
    te, ye, ie, ce = [], [], [], []
    pendingEvents = {}

    # After firing an event, we integrate for root_grace_t without looking for
    # events to prevent finding the same one over and over again due to
    # numerical imprecision
    root_grace_t = (times[-1] - times[0])/1e6
    event_just_fired = False

    # This is the jacobian for use with ddaskr.
    def _ddaskr_jac(t, y, yprime, pd, cj):
        dc_term = net.dres_dc_function(t, y, yprime) 
        dcdot_term = net.dres_dcdot_function(t, y, yprime)
        return dc_term + cj * dcdot_term

    # If calculate_ic is true, then we need to calculate inititial conditions
    # for the DAE.  For daskr to do this it needs to know which variables are
    # algebraic.
    variable_types = None
    if calculate_ic == True:
        variable_types = scipy.zeros(len(IC))
        for jj, vt in enumerate(net._dynamic_var_algebraic):
            variable_types[jj] = vt
   
    exception_raised = None

    while start < times[-1]:

        # check to see if there are any pendingEvents with execution time
        # equal to the start time
        # since some chained events may be added to the pending events list,
        # keep looping until the key has been deleted
        while pendingEvents.has_key(round(start, 12)):

            # Note: the correct behavior of chained events with events that have
            # delays is not clear
            # Note: When multiple events execute simultaneously, we check for
            # chained events after all the simultaneously excecutions occur.
            
            # store a copy of the value of root_func in order to check for
            # chains of events       
            root_before = root_func(net, IC, start).copy()

            # For those events whose execution time == start, execute
            # the event
            for ii in range(len(pendingEvents[round(start, 12)])):
                event = pendingEvents[round(start, 12)].pop()
                IC = net.executeEvent(event, IC, start)
                ypIC = net.get_ddv_dt(IC,start)

            # Update the root state after all listed events have excecuted.
            root_after = root_func(net, IC, start).copy()

            # check for chained events
            dire_cross = scipy.array(root_after) - scipy.array(root_before)
            for ii in range(num_events):
                if dire_cross[ii] > 0 and ii < len(net.events):
                    # an event is firing, so record the time and state
                    te.append(start)
                    ye.append(IC)
                    ie.append(ii)
                    event_just_fired = True
                    event = net.events[ii]
                    logger.debug('Event %s fired at time=%g in network %s.'
                             % (event.id, start, net.id))                              
                    delay = net.fireEvent(net.events[ii], IC, start)
                    if not pendingEvents.has_key(round(start + delay, 12)):
                        pendingEvents[round(start + delay, 12)] = []
                    pendingEvents[round(start + delay, 12)].append(event)
                    root_before[ii] = 0.5

            # if there are no more events to excecute at this time, then
            # delete the entry for this time in pendingEvents
            if len(pendingEvents[round(start, 12)]) == 0:
                del pendingEvents[round(start, 12)]

        # if there are still pending events, set the nextEventTime
        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())
        else:
            nextEventTime = scipy.inf

        # If an event just fired, we integrate for root_grace_t without looking
        # for events, to prevent detecting the same event several times.
        # This section only excecutes if root_grace is True (default).
        if (event_just_fired and root_grace):
            next_requested = scipy.compress(scipy.asarray(times) > start, times)[0]
            integrate_to = min(start + root_grace_t, next_requested,
                               nextEventTime)
            curTimes = [start, integrate_to]

            outputs = integrate_tidbit(net, res_func, _ddaskr_jac, 
                                       root_func=None, 
                                       IC=IC, yp0=ypIC, curTimes=curTimes, 
                                       rtol=rtol, atol=atol, 
                                       fill_traj=fill_traj, 
                                       return_derivs=return_derivs, 
                                       redirect_msgs=redirect_msgs,
                                       init_consistent=False,
                                       var_types=None)



            # only grab the first 4 outputs since we aren't chcking for events
            exception_raised, yout_this, tout_this, youtdt_this = outputs[:4]


            # We don't append the last point, to prevent a needless 'event
            #  looking' duplication of times in the trajectory.
            tout.extend(tout_this[:-1])
            yout = scipy.concatenate((yout, yout_this[:-1]))
            if return_derivs :
                youtdt = scipy.concatenate((youtdt,outputs[3][:-1]))
            event_just_fired = False
            logger.debug('Finished event grace time integration.')

            if exception_raised:
                break

        # If we have pending events, only integrate until the next one.
        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())
            curTimes = scipy.compress((times > start) & (times < nextEventTime),
                                      times)
            curTimes = scipy.concatenate(([start], curTimes, [nextEventTime]))
        else:
            curTimes = scipy.compress(times > start, times)
            curTimes = scipy.concatenate(([start], curTimes))

        outputs = integrate_tidbit(net, res_func, _ddaskr_jac, 
                                   root_func=root_func, 
                                   IC=IC, yp0=ypIC, curTimes=curTimes, 
                                   rtol=rtol, atol=atol, 
                                   fill_traj=fill_traj, 
                                   return_derivs=return_derivs, 
                                   redirect_msgs=redirect_msgs,
                                   init_consistent = calculate_ic,                                         
                                   var_types = variable_types)            

        exception_raised, yout_this, tout_this, youtdt_this,\
              t_root_this, y_root_this, i_root_this = outputs


        yout = scipy.concatenate((yout, yout_this))
        if return_derivs :
            youtdt = scipy.concatenate((youtdt, youtdt_this))
            
        tout.extend(tout_this)

        if exception_raised:
            break

        start, IC = tout[-1], copy.copy(yout[-1])
        ypIC = net.get_ddv_dt(IC,start)

        if return_derivs:
            ypIC = copy.copy(youtdt[-1])
        

        # Check the directions of our root crossings to see whether events
        # actually fired.  DASKR automatically returns the direction of event
        # crossings, so we can read this information from the integration
        # results
        # zip(...) allows us to iterate over multiple sequences in parallel.
        for ii, dire_cross in zip(range(num_events), i_root_this):
            
            # dire_cross is the direction of the event crossing.  If it's
            # 1 that means Ri changed from negative to positive.  In that case
            # we need to record the event.
        
            # We now have each clause of the complex events included in
            # the root_func,  so we need to make sure what we saw was a
            # real event.
            if dire_cross > 0 and ii < len(net.events):
                # an event has fired, so record the state and event number
                te.append(t_root_this)
                ye.append(y_root_this)
                ie.append(ii)
                event_just_fired = True
                event = net.events[ie[-1]]
                logger.debug('Event %s fired at time=%g in network %s.'
                             % (event.id, te[-1], net.id))          
                delay = net.fireEvent(event, yout[-1], te[-1])
                if (pendingEvents.has_key(round(te[-1] + delay, 12)) == False):
                    pendingEvents[round(te[-1] + delay, 12)] = []
                pendingEvents[round(te[-1] + delay, 12)].append(event)

                # handle event subclauses
                # Which clauses of the event fired?
                sub_clauses = scipy.compress(i_root_this >= len(net.events),
                                             i_root_this)                 
                
                event = net.events[ie[-1]]
                if len(sub_clauses) == 0:
                    clause_index = ie[-1]
                elif len(sub_clauses) == 1:
                    clause_index = sub_clauses[0]
                else:
                    raise ValueError('More than one clause in event trigger %s '
                                     'changed simultaneously! Sensitivity '
                                     'integration will fail.' % event.trigger)
                ce.append(clause_index)
            


    # End of while loop for integration.

    if len(yout) and len(tout):
        net.updateVariablesFromDynamicVars(yout[-1], tout[-1])

    if not fill_traj and not exception_raised:
        yout = _reduce_times(yout, tout, times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, times)
        tout = times
    elif reduce_space :
        # make sure we don't miss data times by adding them in.
        # The set usage ensures that we don't end up with duplicates
        filtered_times = [tout[i] for i in range(0,len(tout),reduce_space)]
        filtered_times = sets.Set(filtered_times)
        filtered_times.union_update(times)
        filtered_times = scipy.sort(list(filtered_times))

        yout = _reduce_times(yout, tout, filtered_times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, filtered_times)
        tout = filtered_times
    
    trajectory = Trajectory_mod.Trajectory(net,holds_dt=return_derivs)
    if return_derivs :
        yout = scipy.concatenate((yout,youtdt),axis=1) # join columnwise

    trajectory.appendFromODEINT(tout, yout, holds_dt=return_derivs)
    trajectory.add_event_info((te,ye,ie,ce))
    net.trajectory = trajectory

    # raise exception if exited integration loop prematurely
    if exception_raised:
        logger.warn('Integration ended prematurely in network %s on node %i.'
                    % (net.id, my_rank))
        raise exception_raised

    return trajectory

def integrate_sens_subset(net, times, rtol=1e-6,
                          fill_traj=False, opt_vars=None,
                          return_derivs=False, redirect_msgs=True):
    """
    Integrate the sensitivity equations for a list of optimizable variables.
    """
    rtol, atol = generate_tolerances(net, rtol)
    # We integrate the net once just to know where all our events are
    traj = integrate(net, times, rtol=rtol, fill_traj=fill_traj, 
                     return_events=False, return_derivs=return_derivs,
                     redirect_msgs=redirect_msgs)

    N_dyn_vars = len(net.dynamicVars)
    # Create the array that will hold all our sensitivities
    out_size = (len(traj.get_times()), N_dyn_vars + N_dyn_vars*len(opt_vars))
    all_yout = scipy.zeros(out_size, scipy.float_)

    # Copy the values from the trajectory into this array
    for ii, dyn_var in enumerate(net.dynamicVars.keys()):
        all_yout[:, ii] = traj.get_var_traj(dyn_var)

    # Similarly for the derivative outputs
    if return_derivs:
        all_youtdt = scipy.zeros(out_size, scipy.float_)
        for ii, dyn_var in enumerate(net.dynamicVars.keys()):
            all_youtdt[:, ii] = traj.get_var_traj((dyn_var, 'time'))

    # Now we integrate one optizable variable at a time
    for ii, opt_var in enumerate(opt_vars):
        single_out = integrate_sens_single(net, traj, rtol, opt_var,
                                           return_derivs, redirect_msgs)
        # Copy the result into our main arrays.
        start_col = (ii+1) * N_dyn_vars
        end_col = (ii+2) * N_dyn_vars
        if not return_derivs:
            all_yout[:, start_col:end_col] = single_out
        else:
            all_yout[:, start_col:end_col] = single_out[0]
            all_youtdt[:, start_col:end_col] = single_out[1]

    if not return_derivs:
        return traj.get_times(), all_yout
    else:
        return traj.get_times(), all_yout, all_youtdt

def integrate_sens_single(net, traj, rtol, opt_var, return_derivs,
                          redirect_msgs):
    """
    Integrate the sensitivity equations for a single optimization variable.
    """
    #
    # This now uses the ddaskr related-functions.
    #
    opt_var_index = net.optimizableVars.index_by_key(opt_var)
    N_dyn_vars = len(net.dynamicVars)

    # Calculate tolerances for our sensitivity integration
    rtol, atol = generate_tolerances(net, rtol)
    # We set the same rtol for sensitivity variables as for the normal
    # variables.
    rtol_for_sens = scipy.concatenate((rtol, rtol))
    # Absolute tolerances depend on the typical value of the optimizable
    # variable
    atol_for_sens = atol/net.get_var_typical_val(opt_var)
    atol_for_sens = scipy.concatenate((atol, atol_for_sens))

    # The passed-in normal trajectory tells us what times we need to return,
    #  and what events will happen.
    times = traj.get_times()
    te, ye, ie, ce = traj.event_info

    current_time = times[0]

    # Initial conditions
    IC = scipy.zeros(2*N_dyn_vars, scipy.float_)
    for dvInd, (id, var) in enumerate(net.dynamicVars.items()):
        # For normal integration
        IC[dvInd] = traj.get_var_val(id, current_time)
        # For sensitivity integration
        if isinstance(var.initialValue, str):
            DwrtOV = net.takeDerivative(var.initialValue, opt_var)
            IC[N_dyn_vars + dvInd] = net.evaluate_expr(DwrtOV, 
                                                       time=current_time)
    # We'll let ddaskr figure out the initial condition on yprime
    ypIC = 0*IC + 1
    # The var_types for the sensitivity variables are identical to those
    # for the normal variables.
    var_types = scipy.concatenate((net._dynamic_var_algebraic,
                                   net._dynamic_var_algebraic))

    tout = []
    sens_out = scipy.zeros((0, N_dyn_vars), scipy.float_)
    sensdt_out = scipy.zeros((0, N_dyn_vars), scipy.float_)

    pendingEvents = {}
    while current_time < times[-1]:
        if round(current_time, 12) in pendingEvents.keys():
            IC = net.executeEventAndUpdateDdv_Dov(event, IC, current_time, 
                                                  [opt_var])
            del pendingEvents[round(current_time, 12)]
        # Figure out how long we need to integrate until
        if pendingEvents:
            integrate_until = min(pendingEvents.keys())
        else:
            integrate_until = scipy.inf
        if te:
            integrate_until = min(integrate_until, min(te))

        # The the list of times in this part of the integration
        int_times = scipy.compress((scipy.asarray(times) > current_time)
                                  & (scipy.asarray(times) < integrate_until),
                                  times)
        int_times = scipy.concatenate(([current_time], int_times))
        if integrate_until < times[-1]:
            int_times = scipy.concatenate((int_times, [integrate_until]))

        # We let daskr figure out ypIC
        if net.integrateWithLogs:
            IC[:N_dyn_vars] = scipy.log(IC[:N_dyn_vars])
        try:
            int_outputs = daeint(_alg_sens, int_times,
                                 IC, ypIC, 
                                 rtol_for_sens, atol_for_sens,
                                 args = (net, opt_var_index),
                                 max_steps = 1e4,
                                 init_consistent=1,
                                 var_types = var_types,
                                 redir_output = redirect_msgs)
        except Utility.SloppyCellException, X:
            logger.warn('Sensitivity integration failed for network %s on '
                        'node %i during optimizable variable %s.'
                        % (net.id, my_rank, opt_var))
            raise

        # Copy out the sensitivity values
        yout_this = int_outputs[0]
        youtdt_this = int_outputs[2]
        sens_out = scipy.concatenate((sens_out, 
                                      yout_this[:,N_dyn_vars:].copy()))
        if return_derivs:
            sensdt_out = scipy.concatenate((sensdt_out,
                                            youtdt_this[:,N_dyn_vars:].copy()))
        tout.extend(int_times)

        current_time, IC = tout[-1], yout_this[-1].copy()
        
        # All that matters is the IC
        if net.integrateWithLogs:
            IC[:N_dyn_vars] = scipy.exp(IC[:N_dyn_vars])

        # If an event should be fired
        if te and te[0] == tout[-1]:
            event = net.events[ie[0]]
            logger.debug('Event %s fired at time=%g in network %s for sens.'
                         % (event.id, te[0], net.get_id()))
            delay = net.fireEventAndUpdateDdv_Dov(event, IC, te[0], 
                                                  net.event_clauses[ce[0]],
                                                  [opt_var])
            pendingEvents[round(te[0] + delay, 12)] = event
            # Remove this event from the lists of events
            te = te[1:]
            ie = ie[1:]
            ce = ce[1:]

    sens_out = _reduce_times(sens_out, tout, times)
    if not return_derivs:
        return sens_out
    else:
        sensdt_out = _reduce_times(sensdt_out, tout, times)
        return sens_out, sensdt_out

def _parse_sens_result(result, net, opt_vars, yout, youtdt=None):
    """
    Utility function for parsing the return from integrate_sens_subset
    """
    n_dyn, n_opt = len(net.dynamicVars), len(net.optimizableVars)
    for work_index, opt_var in enumerate(opt_vars):
        net_index = net.optimizableVars.indexByKey(opt_var)
        net_start = (net_index+1)*n_dyn
        net_end = (net_index+2)*n_dyn

        work_start = (work_index+1)*n_dyn
        work_end = (work_index+2)*n_dyn

        yout[:, net_start:net_end] = result[1][:, work_start:work_end]
        if youtdt is not None:
            youtdt[:, net_start: net_end] =\
                    result[2][:, work_start:work_end]

def integrate_sensitivity(net, times, params=None, rtol=1e-6, 
                          fill_traj=False, return_derivs=False,
                          redirect_msgs=True):
    logger.debug('Entering integrate_sens on node %i' % my_rank)

    times = scipy.array(times)
    net.compile()
    if times[0] == 0:
        net.resetDynamicVariables()

    if params is not None:
        net.update_optimizable_vars(params)

    # Assigned variables to process on each node.
    vars_assigned = dict([(node, net.optimizableVars.keys()[node::num_procs]) 
                            for node in range(num_procs)])

    # Send out jobs for workers
    for worker in range(1, num_procs):
        logger.debug('Sending to worker %i: %s' % (worker, 
                                                   str(vars_assigned[worker])))
        command = 'Dynamics.integrate_sens_subset(net, times,'\
                'rtol, fill_traj, opt_vars, return_derivs, redirect_msgs=redir)'
        args = {'net':net, 'times':times, 'rtol':rtol, 'fill_traj':fill_traj,
                'opt_vars':vars_assigned[worker], 'return_derivs':return_derivs,
                'redir': redirect_msgs}
        pypar.send((command, args), worker)

    logger.debug('Master doing vars %s' % str(vars_assigned[0]))
    try:
        result = integrate_sens_subset(net, times, rtol, fill_traj, 
                                       vars_assigned[0], return_derivs,
                                       redirect_msgs=redirect_msgs)
    except Utility.SloppyCellException:
        # If the master encounters an exception, we must still wait and get 
        #  replies from all the workers (even if we do nothing with them) so 
        #  that communication stays synchronized.
        for worker in range(1, num_procs):
            pypar.receive(worker)
        raise

    # Begin pulling results together...
    tout = result[0]
    # Build the yout array
    n_dyn, n_opt = len(net.dynamicVars), len(net.optimizableVars)
    yout = scipy.zeros((len(tout), n_dyn * (n_opt+ 1)), scipy.float_)

    # We use the master's result for the non-sensitivity values
    yout[:, :n_dyn] = result[1][:, :n_dyn]
    if return_derivs:
        youtdt = scipy.zeros((len(tout), n_dyn * (n_opt+ 1)), scipy.float_)
        youtdt[:, :n_dyn] = result[2][:, :n_dyn]
    else:
        youtdt = None

    # Copy the sensitivity results into yout and (if necessary) youtdt
    _parse_sens_result(result, net, vars_assigned[0], yout, youtdt)

    # Now when we listen to the worker's replies, we store any exception they
    #  return in exception_raised. We'll reraise that exception after getting
    #  replies from all the workers.
    exception_raised = None
    for worker in range(1, num_procs):
        logger.debug('Receiving result from worker %i.' % worker)
        result = pypar.receive(worker)
        if isinstance(result, Utility.SloppyCellException):
            exception_raised = result
            continue
        _parse_sens_result(result, net, vars_assigned[worker], yout, youtdt)

    if exception_raised:
        raise exception_raised

    ddv_dpTrajectory = Trajectory_mod.Trajectory(net, is_sens=True, 
                                                 holds_dt=return_derivs)
    if return_derivs:
        yout = scipy.concatenate((yout, youtdt), axis=1)
    ddv_dpTrajectory.appendSensFromODEINT(tout, yout, holds_dt = return_derivs)

    net.trajectory = ddv_dpTrajectory

    return ddv_dpTrajectory

def dyn_var_fixed_point(net, dv0=None, with_logs=True, xtol=1e-6):
    """
    Return the dynamic variables values at the closest fixed point of the net.

    dv0  Initial guess for the fixed point. If not given, the current state
         of the net is used.
    with_logs   If True, the calculation is done in terms of logs of variables,
                so that they cannot be negative.
    """
    net.compile()

    if dv0 is None:
        dv0 = scipy.array(net.getDynamicVarValues())
        # We take the absolute value of dv0 to avoid problems from small 
        #  numerical noise negative values.
        if scipy.any(dv0 < 0):
            logger.warning('Negative values in initial guess for fixed point. '
                           'Taking absolute values. The most negative value '
                           'was %g.' % min(dv0))
            dv0 = scipy.absolute(dv0)

    if with_logs:
        func = lambda logDV: net.get_ddv_dt(scipy.exp(logDV), 0)
        def fprime(logDV):
            fp = net.get_d2dv_ddvdt(scipy.exp(logDV), 0)
            # Our derivatives run down the columns, so we need to 
            #  take the tranpose twice to multiply down them.
            return scipy.transpose(scipy.transpose(fp) * scipy.exp(logDV))
        x0 = scipy.log(dv0)
        # To transform sigma_x to sigma_log_x, we divide by x. We can set
        #  our sigma_log_x use to be the mean of what our xtol would yield
        #  for each chemical.
        xtol = scipy.mean(xtol/dv0)
    else:
        func = lambda dv: net.get_ddv_dt(dv, 0)
        fprime = lambda dv: net.get_d2dv_ddvdt(dv, 0)
        x0 = dv0

    try:
        dvFixed, infodict, ier, mesg =\
                scipy.optimize.fsolve(func, x0=x0, fprime=fprime,
                                      col_deriv=True, full_output=True,
                                      xtol=xtol)
    except (scipy.optimize.minpack.error, ArithmeticError), X:
        raise FixedPointException(('Failure in fsolve.', X))

    if ier != 1:
        raise FixedPointException(mesg, infodict)

    if with_logs:
        return scipy.exp(dvFixed)
    else:
        return dvFixed

def _alg_sens(time, y, ydot, net, opt_ii):
    # This is the integrand for the sensitivity integration
    nDV = len(net.dynamicVars)

    if not net.integrateWithLogs:
        c = y[:nDV]
        cdot = ydot[:nDV]
    else:
        c = scipy.exp(y[:nDV])
        # Force our c to be non-zero, since that's the whole point of the
        # logarithmic integration.
        c = scipy.maximum(c, scipy.misc.limits.double_tiny)
        cdot = ydot[:nDV] * c
    dc_dp = y[nDV:]
    dcdot_dp = ydot[nDV:]

    # These are the residuals for the normal integration
    res = net.res_function(time, c, cdot, ires=0)

    # The total derivative of res with respect to parameter opt_ii is
    # Dres_Dp = \partial res/\partial dp + \partial res/\partial c * dc/dp +
    #    \partial res/\partial cdot * dcdot/dp
    dres_dp = net.dres_dp_function(time, c, cdot, 
                                   indices=[opt_ii])[:, opt_ii]
    dres_dc = net.dres_dc_function(time, c, cdot)
    dres_dcdot = net.dres_dcdot_function(time, c, cdot)

    Dres_Dp = dres_dp + scipy.dot(dres_dc, dc_dp) +\
            scipy.dot(dres_dcdot, dcdot_dp)

    # Glue them together for the return
    result = scipy.concatenate((res, Dres_Dp))
    return result


def _reduce_times(yout, tout, times):
    jj = 0
    maxjj = len(tout)
    for ii, twanted in enumerate(times):
        while (tout[jj] != twanted) & (jj<maxjj-1):
            jj += 1
        yout[ii] = yout[jj]

    yout = yout[:len(times)]

    return yout
