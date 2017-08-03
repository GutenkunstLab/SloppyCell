"""
Methods for evaluating the dynamics of Network.
"""
__docformat__ = "restructuredtext en"

import copy
import sets
import sys 
import scipy
import scipy.optimize

import logging
logger = logging.getLogger('ReactionNetworks.Dynamics')

import SloppyCell.daskr as daskr
daeint = daskr.daeint
import Trajectory_mod

import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks.Components import event_info
from SloppyCell import HAVE_PYPAR, my_rank, my_host, num_procs
if HAVE_PYPAR:
    import pypar

_double_epsilon_ = scipy.finfo(scipy.float_).eps
_double_tiny_ = scipy.finfo(scipy.float_).tiny

# This specifies the maximum number of steps daeint will take between timepoints
MAX_STEPS = 1e5

# The maximum number of timepoints daeint is allowed to return.
MAX_TIMEPOINTS = 1e5

# global_rtol is the default relative tolerance to use if none is specied in
# calls.
global_rtol = 1e-6
# global_atol overrides the default of rtol*typical_value. It may be useful if
# we're having problems with large values going negative. (Although then it
# might be better to just use logarithmic integration.)
global_atol = None
global_hmax = None

return_derivs = False # do we want time derivatives of all trajectories 
                      # returned?
reduce_space = 0 # we may want to not return every timepoint but only 0, 
                 # reduce_space, 2*reduce_space etc.

failed_args = []
failed_kwargs = {}

class IntegrationException(Utility.SloppyCellException):
    pass

class FixedPointException(Utility.SloppyCellException):
    pass

def integrate_tidbit(net, res_func, Dfun, root_func, IC, yp0, curTimes, 
                     rtol, atol, fill_traj, return_derivs,
                     redirect_msgs, calculate_ic, var_types):
    N_dyn_var = len(net.dynamicVars)

    if net.integrateWithLogs:
        yp0 = yp0/IC
        IC = scipy.log(IC)

        res_func = net.res_function_logdv
        Dfun = None
        if root_func is not None:
            root_func = net.root_func_logdv

        atol = rtol
        rtol = [0]*len(rtol)

    if root_func is None:
        nrt = 0
    else:
        nrt = net.len_root_func

    int_args = {'res': res_func,
                't': curTimes,
                'y0': IC,
                'yp0': yp0,

                'jac': Dfun,
                'nrt': nrt,
                'rt': root_func,

                'rpar': net.constantVarValues,

                'rtol': rtol,
                'atol': atol,
                'max_steps': MAX_STEPS,

                'intermediate_output': fill_traj,
                'redir_output': redirect_msgs,

                'calculate_ic' : calculate_ic,
                'var_types' : var_types,
                'hmax' : global_hmax,
        
                'max_timepoints': MAX_TIMEPOINTS
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
        ydt_this = ydt_this * y_this
        y_root_this = scipy.exp(y_root_this)

    return exception_raised, y_this, t_this, ydt_this,\
            t_root_this, y_root_this, i_root_this

def generate_tolerances(net, rtol, atol=None):
    if rtol is None:
        rtol = global_rtol
    if scipy.isscalar(rtol):
        rtol = scipy.ones(len(net.dynamicVars)) * rtol

    # We set atol to be a minimum of global_atol to avoid variables with large
    # typical values going negative.
    if (scipy.isscalar(atol) or atol==None):
        typ_vals = [abs(net.get_var_typical_val(id))
                    for id in net.dynamicVars.keys()]
        atol = rtol * scipy.asarray(typ_vals)
        if global_atol:
            atol = scipy.minimum(atol, global_atol)
    return rtol, atol


def integrate(net, times, rtol=None, atol=None, params=None, fill_traj=True,
              return_events=False, return_derivs=False,
              redirect_msgs=True, calculate_ic = False,
              include_extra_event_info = False, use_constraints=True):
    """
    Integrate a Network, returning a Trajectory.

    net            The Network to integrate.
    times          A sequence of times to include in the integration output.
    params         Parameters for the integration. If params=None, the current
                   parameters of net are used.
    rtol           Relative error tolerance for this integration.  
                   If rtol = None then relative tolerance is set by 
                   Dynamics.global_rtol.
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
    include_extra_event_info    If True, the returned trajectory will have more
                    detailed event information, inclusing pre- and post- execution
                    state, and assigned variable informtaion.
    use_constraints    If True, and if the network has constraints, will raise
                       an exception if a constraint's math function becomes 
                       False at any time.
    """
    logger.debug('Integrating network %s.' % net.get_id())

    if params is not None:
        net.update_optimizable_vars(params)
    constants = net.constantVarValues
    # On some systems, f2py'd functions fail if len(constants) == 0.
    if len(constants) == 0:
        constants = [0]
    # If you ask for time = 0, we'll assume you want dynamic variable values
    # reset.
    times = scipy.asarray(times)
    if times[0] == 0:
        net.resetDynamicVariables()
    net.compile()

    if (rtol is None or atol is None):
        rtol, atol = generate_tolerances(net, rtol, atol)

    # get the initial state, time, and derivative.
    start = times[0]
    IC = net.getDynamicVarValues()
    # We start with ypIC equal to typ_vals so that we get the correct scale.
    typ_vals = [net.get_var_typical_val(id) for id in net.dynamicVars.keys()]
    ypIC = scipy.array(typ_vals)
    # This calculates the yprime ICs and also ensures that values for our
    #  algebraic variables are all consistent.
    IC, ypIC = find_ics(IC, ypIC, start,
                        net._dynamic_var_algebraic, rtol, atol,
                        net.constantVarValues, net, 
                        redirect_msgs=redirect_msgs)

    """

1) After the initial conditions have been calculated, use
net.updateVariablesFromDynamicVars to make sure the Network object has
the proper variable values.
2) Then loop through the constraint events.
    a) For each constraint, result = net.evaluate_expr(constraint.trigger).
    b) If result is a violation (I don't know whether that means True
or False), raise the exception.
    """
    # check that constraints are not violated at t=0
    if use_constraints == True:
        net.updateVariablesFromDynamicVars(values=IC, time=start)
        for con_id, constraint in net.constraints.items():
            result = net.evaluate_expr(constraint.trigger)
            if result == False:
                raise Utility.ConstraintViolatedException(start,
                                                 constraint.trigger,
                                                 constraint.message)
                
            
    
    # start variables for output storage 
    yout = scipy.zeros((0, len(IC)), scipy.float_)
    youtdt = scipy.zeros((0, len(IC)), scipy.float_)
    tout = []

    # For some reason, the F2Py wrapper for ddaskr seems to fail if this
    #  isn't wrapped up in a lambda...
    res_func = net.res_function
    root_func = net.root_func

    # events_occured will hold event_info objects that describe the events
    #  in detail.
    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events. Here mostly for historical reasons.
    events_occurred = []
    te, ye, ie = [], [], []
    pendingEvents = {}

    event_just_fired = False
    event_just_executed = False

    # This is the jacobian for use with ddaskr.
    try:
        _ddaskr_jac = net.ddaskr_jac
    except AttributeError:
        _ddaskr_jac = None

    exception_raised = None
    while start < times[-1]:
        # check to see if there are any pendingEvents with execution time
        # equal to the start time
        # since some chained events may be added to the pending events list,
        # keep looping until the key has been deleted
        if pendingEvents and min(pendingEvents.keys()) < start:
            raise ValueError('Missed an event!')
        event_buffer = 0
        while pendingEvents.has_key(start):
            execution_time = start
            # We need to backtrack to deal with this event...
            event_list = pendingEvents[execution_time]

            # Note: the correct behavior of chained events with events that have
            # delays is not clear
            # Note: When multiple events execute simultaneously, we check for
            # chained events after all the simultaneously excecutions occur.
            
            # store a copy of the value of root_func in order to check for
            # chains of events       
            root_before = root_func(start, IC, ypIC, constants)

            # For those events whose execution time == start, execute
            # the event
            # We also record the event that started this chain of events firing.
            # This is necessary for sensitivity integration.
            chain_starter = event_list[-1]
            while len(event_list) > 0:
                holder = event_list.pop()

                # check whether the current event was a constraint
                # technically this should be caught in fired_events, so the
                # if statement below is deprecated and can be removed.
                if holder.event_index >= len(net.events) and use_constraints == True:
                    # a constraint was violated
                    con_id = holder.event_id
                    constraint = net.constraints.get(con_id)
                    raise Utility.ConstraintViolatedException(start,
                                                     constraint.trigger,
                                                     constraint.message)
                holder.time_exec = start
                holder.y_pre_exec = copy.copy(IC)
                holder.yp_pre_exec = copy.copy(ypIC)
                IC = net.executeEvent(event=holder.event, 
                                      time_fired=start,
                                      y_fired=holder.y_fired, 
                                      time_current=start,
                                      y_current=IC)
                # Now we update the initial conditions to reflect any changes
                #  made in event execution. This might include adjustments
                #  to algebraic variables and to yprime.
                IC, ypIC = find_ics(IC, ypIC, start,
                                    net._dynamic_var_algebraic, rtol, atol,
                                    net.constantVarValues, net,
                                    redirect_msgs=redirect_msgs)
                holder.y_post_exec = copy.copy(IC)
                holder.yp_post_exec = copy.copy(ypIC)
                event_buffer = max(event_buffer, holder.event.buffer)

            if event_buffer:
                # Ensure that we don't event_buffer past any
                # requested timepoints.
                nextEventTime = scipy.inf
                pendingEventTimes = pendingEvents.keys()
                pendingEventTimes.remove(execution_time)
                if pendingEventTimes:
                    nextEventTime = min(pendingEventTimes)
                else:
                    nextEventTime = scipy.inf
                nextEventTime = min([nextEventTime, start+event_buffer])
                if nextEventTime < times[-1]:
                    curTimes = scipy.compress((times > start) 
                                              & (times < nextEventTime), times)
                    curTimes = scipy.concatenate(([start], curTimes, 
                                                  [nextEventTime]))
                else:
                    curTimes = scipy.compress(times > start, times)
                    curTimes = scipy.concatenate(([start], curTimes))
                outputs = integrate_tidbit(net, res_func, _ddaskr_jac, 
                                           root_func=None, 
                                           IC=IC, yp0=ypIC, 
                                           curTimes=curTimes,
                                           rtol=rtol, atol=atol, 
                                           fill_traj=fill_traj, 
                                           return_derivs=True, 
                                           redirect_msgs=redirect_msgs,
                                           calculate_ic = False,
                                           var_types=net._dynamic_var_algebraic)

                exception_raised, yout_this, tout_this, youtdt_this,\
                      t_root_this, y_root_this, i_root_this = outputs

                yout = scipy.concatenate((yout, yout_this))
                youtdt = scipy.concatenate((youtdt, youtdt_this))
                tout.extend(tout_this)

                if exception_raised:
                    break

                start = tout[-1]
                IC = copy.copy(yout[-1])
                ypIC = copy.copy(youtdt_this[-1])
            # We need to break several times to get to the proper level to deal
            # with an exception.
            if exception_raised:
                break

            # Update the root state after all listed events have excecuted.
            root_after = root_func(start, IC, ypIC, constants)

            # Check for chained events/constraints
            crossing_dirs = root_after - root_before
            event_just_fired = fired_events(net, start, IC, ypIC, 
                                            crossing_dirs,
                                            events_occurred, pendingEvents,
                                            chained_off_of = chain_starter,
                                            use_constraints=use_constraints)
            event_just_executed = True

            # If there are no more events to excecute at this time, then
            #  delete the entry for this time in pendingEvents.
            # Otherwise we go back to the top of this loop, to excute the
            #  next event, check for firing events, etc.
            if len(pendingEvents[execution_time]) == 0:
                del pendingEvents[execution_time]
        # We need to break several times to get to the proper level to deal
        # with an exception.
        if exception_raised:
            break

        # We remove 'double timepoints' that we would otherwise have when
        #  an event fired, but didn't execute due to a delay.
        if not event_just_executed and len(tout) > 0:
            tout = tout[:-1]
            yout = yout[:-1]
            youtdt = youtdt[:-1]
        event_just_executed = False

        # If there are still pending events, set the nextEventTime
        nextEventTime = scipy.inf
        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())

        # If we have pending events, only integrate until the next one.
        if nextEventTime < times[-1]:
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
                                   return_derivs=True, 
                                   redirect_msgs=redirect_msgs,
                                   calculate_ic = False,
                                   var_types = net._dynamic_var_algebraic)

        exception_raised, yout_this, tout_this, youtdt_this,\
              t_root_this, y_root_this, i_root_this = outputs

        yout = scipy.concatenate((yout, yout_this))
        youtdt = scipy.concatenate((youtdt, youtdt_this))
        tout.extend(tout_this)

        if exception_raised:
            break

        start = tout[-1]
        IC = copy.copy(yout[-1])
        ypIC = copy.copy(youtdt_this[-1])

        # Check for events firing.
        # If one of the fired events is a constraint, then we'll catch
        # it at the top of the while loop for integration
        event_just_fired = fired_events(net, start, IC, ypIC, 
                                                i_root_this, 
                                                events_occurred, pendingEvents,
                                                use_constraints=use_constraints)

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
    
    trajectory = Trajectory_mod.Trajectory(net, holds_dt=return_derivs,
                                           const_vals=net.constantVarValues)
    if return_derivs :
        yout = scipy.concatenate((yout,youtdt),axis=1) # join columnwise

    trajectory.appendFromODEINT(tout, yout, holds_dt=return_derivs)
    # Fill in the historical te, ye, ie lists
    te, ye, ie, ye_post = [],[],[],[]
    for holder in events_occurred:
        # Here we only record those events that executed. It is possible to fire
        # during a simulation (so be in events_occurred), but not execute due
        # to a delay or to being at the very end of the integration.
        if hasattr(holder, 'time_exec'):
            te.append(holder.time_exec)
            ye.append(holder.y_pre_exec)
            ye_post.append(holder.y_post_exec)
            ie.append(holder.event_index)

    # We only include y_post_exec and assignned vars if requested
    if include_extra_event_info == False:
        trajectory.add_event_info(net, (te,ye,ie), tout[-1],
                              include_extra_event_info = include_extra_event_info)
    elif include_extra_event_info == True:
        trajectory.add_event_info(net, (te,ye,ye_post,ie), tout[-1],
                              include_extra_event_info = include_extra_event_info)


    trajectory.events_occurred = events_occurred
    net.trajectory = trajectory

    # raise exception if exited integration loop prematurely
    if exception_raised:
        logger.warn('Integration ended prematurely in network %s on node %i.'
                    % (net.id, my_rank))
        raise exception_raised

    return trajectory

def fired_events(net, time, y, yp, crossing_dirs, 
                 events_occurred, pendingEvents,
                 chained_off_of=None, use_constraints=False):
    # checking for both normal events and constraint events
    num_events = len(net.events)+len(net.constraints)
    event_just_fired = False
    # Check the directions of our root crossings to see whether events
    # actually fired.  DASKR automatically returns the direction of event
    # crossings, so we can read this information from the integration
    # results.
    # Note that we might have more than one event firing at the same time.

    for event_index, dire_cross in zip(range(num_events), crossing_dirs):
        # dire_cross is the direction of the event crossing.  If it's
        # 1 that means Ri changed from negative to positive.  In that case
        # we need to record the event. Note that num_events runs over the
        # real events, not the sub-clauses which are also in root_func.
        if dire_cross > 0 and event_index < len(net.events):
            event_just_fired = True
            if event_index < len(net.events):
                event = net.events[event_index]
                
            logger.debug('Event %s fired at time=%g in network %s.'
                         % (event.id, time, net.id))

            # Which clauses of the event fired?
            # These are the crossing_dirs's corresponding to the subclauses
            sub_clause_dirs = crossing_dirs[num_events:]
            # We need to add back in num_events for counting purposes
            sub_clauses_fired = scipy.nonzero(sub_clause_dirs > 0)[0]\
                    + num_events
            # Now we check which, if any, of these sub_clauses actually
            #  belong to this event. This is only necessary because events
            #  may fire at the same time.
            sub_clauses = []
            for fired_index in sub_clauses_fired:
                if fired_index in event.sub_clause_indices:
                    sub_clauses.append(fired_index)
                
            # If no sub-clauses fired, then it was just a simple event.
            if len(sub_clauses) == 0:
                clause_index = event_index
            # But if one of the subclauses fired, we want to record that.
            elif len(sub_clauses) == 1:
                clause_index = sub_clauses.pop()
            else:
                # We might have a problem.
                raise ValueError('More than one clause in event trigger %s '
                                 'changed simultaneously! Sensitivity '
                                 'integration will fail.' % event.trigger)
            delay = net.fireEvent(event, y, time)
            execution_time = time + delay

            # We use this object to record which events have fired
            holder = event_info()
            holder.event = event
            holder.event_index = event_index
            holder.event_id = event.id
            holder.clause_index = clause_index
            holder.time_fired = time
            holder.y_fired = copy.copy(y)
            holder.yp_fired = copy.copy(yp)
            holder.delay = delay
            holder.chained_off_of = chained_off_of
            events_occurred.append(holder)

            # Add this event to the list of events that are supposed to
            # execute.
            if not pendingEvents.has_key(execution_time):
                pendingEvents[execution_time] = []
            pendingEvents[execution_time].append(holder)

        if dire_cross < 0 and event_index >= len(net.events) \
           and event_index <= len(net.events) + len(net.constraints):
            event_just_fired = True
            event = net.constraints[event_index-len(net.events)]
            con_id = event.id
            constraint = net.constraints.get(con_id)
            if use_constraints == True:
                raise Utility.ConstraintViolatedException(time,
                                                 constraint.trigger,
                                                 constraint.message)


    return event_just_fired

def integrate_sens_subset(net, times, rtol=None,
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
            all_yout[:, start_col:end_col] = single_out[0]
        else:
            all_yout[:, start_col:end_col] = single_out[0]
            all_youtdt[:, start_col:end_col] = single_out[1]
        # Concatenate the various 'events_occurred'
        for eii,e in enumerate(traj.events_occurred):
            if not hasattr(e, 'time_exec'):
                # This is an event that fired but never executed. We can ignore
                # it since it doesn't affect the dynamics.
                continue
            if not hasattr(e, 'ysens_fired'):
                e.ysens_fired = single_out[-1][eii].ysens_fired
                e.ysens_post_exec = single_out[-1][eii].ysens_post_exec
                e.ysens_pre_exec = single_out[-1][eii].ysens_pre_exec
            else:
                e.ysens_fired = scipy.concatenate((e.ysens_fired, single_out[-1][eii].ysens_fired[N_dyn_vars:]))
                e.ysens_post_exec = scipy.concatenate((e.ysens_post_exec, single_out[-1][eii].ysens_post_exec[N_dyn_vars:]))
                e.ysens_pre_exec = scipy.concatenate((e.ysens_pre_exec, single_out[-1][eii].ysens_pre_exec[N_dyn_vars:]))

    if not return_derivs:
        return traj.get_times(), all_yout, traj.event_info, traj.events_occurred
    else:
        return traj.get_times(), all_yout, all_youtdt, traj.event_info, traj.events_occurred

def integrate_sens_single(net, traj, rtol, opt_var, return_derivs,
                          redirect_msgs):
    """
    Integrate the sensitivity equations for a single optimization variable.
    """
    logger.debug('Integrating sensitivity for %s.' % opt_var)
    opt_var_index = net.optimizableVars.index_by_key(opt_var)
    N_dyn_vars = len(net.dynamicVars)

    # Calculate tolerances for our sensitivity integration
    rtol, atol = generate_tolerances(net, rtol)
    if net.integrateWithLogs:
        atol = rtol
        rtol = [0] * len(net.dynamicVars)
    # We set the same rtol for sensitivity variables as for the normal
    # variables.
    rtol_for_sens = scipy.concatenate((rtol, rtol))
    # Absolute tolerances depend on the typical value of the optimizable
    # variable
    atol_for_sens = atol/abs(net.get_var_typical_val(opt_var))
    atol_for_sens = scipy.concatenate((atol, atol_for_sens))

    # The passed-in normal trajectory tells us what times we need to return,
    #  and what events will happen.
    times = traj.get_times()
    times = scipy.asarray(times)
    events_occurred = copy.deepcopy(traj.events_occurred)

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

    # The index of the variable to calculate sensitivity wrt is passed in
    # the last element of rpar. It is obtained by an (int) cast, so we add
    # a small amount to make sure the cast doesn't round down stupidly.
    rpar = scipy.concatenate((net.constantVarValues, [opt_var_index+0.1]))

    # Our intial guess for ypIC will be based on the typical values
    typ_vals = [net.get_var_typical_val(id) for id in net.dynamicVars.keys()]*2
    ypIC = scipy.array(typ_vals, dtype=float)
    ypIC[N_dyn_vars:] /= net.get_var_ic(opt_var)

    tout = []
    sens_out = scipy.zeros((0, N_dyn_vars), scipy.float_)
    sensdt_out = scipy.zeros((0, N_dyn_vars), scipy.float_)

    fired_events = {}
    executed_events = {}
    for holder in events_occurred:
        # We only care about events that actually executed.
        if hasattr(holder, 'time_exec'):
            firing_time = holder.time_fired
            fired_events.setdefault(firing_time, [])
            fired_events[firing_time].append(holder)

            execution_time = holder.time_exec
            executed_events.setdefault(execution_time, [])
            executed_events[execution_time].append(holder)

    event_just_executed = False
    while current_time < times[-1]:
        logger.debug('sens_int current_time: %f' % current_time)
        if executed_events and min(executed_events.keys()) < current_time:
            raise ValueError('Missed an event execution in sensitivity '
                             'integration!')
        if executed_events.has_key(current_time):
            event_list = executed_events[current_time]
            for holder in event_list:
                logger.debug('Executing event at time %f.' % current_time)
                IC = net.executeEventAndUpdateSens(holder, IC, opt_var)
            del executed_events[current_time]
            event_just_executed = True

        # We remove 'double timepoints' that we would otherwise have when
        #  an event fired, but didn't execute due to a delay.
        if not event_just_executed and len(tout) > 0:
            tout = tout[:-1]
            sens_out = sens_out[:-1]
            sensdt_out = sensdt_out[:-1]
        event_just_executed = False

        integrate_until = times[-1]
        if fired_events:
            integrate_until = min(integrate_until, min(fired_events.keys()))
        if executed_events:
            integrate_until = min(integrate_until, min(executed_events.keys()))

        # The the list of times in this part of the integration
        if integrate_until < times[-1]:
            int_times = scipy.compress((times > current_time)
                                      & (times < integrate_until),
                                      times)
            int_times = scipy.concatenate(([current_time], int_times, 
                                          [integrate_until]))
        else:
            int_times = scipy.compress(times > current_time, times)
            int_times = scipy.concatenate(([current_time], int_times))

        ypIC = find_ypic_sens(IC, ypIC, current_time, 
                              net._dynamic_var_algebraic,
                              rtol, atol_for_sens, rpar, net, opt_var)

        sens_rhs = net.sens_rhs
        if net.integrateWithLogs:
            ypIC[:N_dyn_vars] *= IC[:N_dyn_vars]
            IC[:N_dyn_vars] = scipy.log(IC[:N_dyn_vars])
            sens_rhs = net.sens_rhs_logdv

        # Note that IC is not necessarily correct for the sensitivities of
        # algebraic variables. We'll let ddaskr calculate them using
        # calculate_ic = True.
        # They will be correct in the returned trajectory. For now, we don't
        # use d/dt of sensitivities, so it doesn't matter if those are wrong
        # at single points in the trajectory.
        var_types = scipy.concatenate((net._dynamic_var_algebraic,
                                       net._dynamic_var_algebraic))
        try:
            int_outputs = daeint(sens_rhs, int_times,
                                 IC, ypIC, 
                                 rtol_for_sens, atol_for_sens,
                                 rpar = rpar,
                                 max_steps = MAX_STEPS,
                                 var_types = var_types,
                                 calculate_ic = True,
                                 redir_output = redirect_msgs,
                                 hmax = global_hmax,
                                 max_timepoints = MAX_TIMEPOINTS)
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
        sensdt_out = scipy.concatenate((sensdt_out,
                                        youtdt_this[:,N_dyn_vars:].copy()))
        tout.extend(int_times)

        current_time, IC = tout[-1], yout_this[-1].copy()
        ypIC = youtdt_this[-1]
        if net.integrateWithLogs:
            IC[:N_dyn_vars] = scipy.exp(IC[:N_dyn_vars])
            ypIC[:N_dyn_vars] *= IC[:N_dyn_vars]

        if fired_events and min(fired_events.keys()) < current_time:
            raise ValueError('Missed event firing in sensitivity integration!')
        if fired_events.has_key(current_time):
            logger.debug('Firing event at time %f.' % current_time)
            event_list = fired_events[current_time]
            while event_list:
                holder = event_list.pop()
                holder.ysens_fired = copy.copy(IC)
            del fired_events[current_time]

    sens_out = _reduce_times(sens_out, tout, times)
    if not return_derivs:
        return sens_out, events_occurred
    else:
        sensdt_out = _reduce_times(sensdt_out, tout, times)
        return sens_out, sensdt_out, events_occurred

def _parse_sens_result(result, net, opt_vars, yout, youtdt=None, events_occurred=None):
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
    if events_occurred is not None:
        for eii,e in enumerate(events_occurred):
            e.ysens_fired = scipy.concatenate((e.ysens_fired, result[-1][eii].ysens_fired[n_dyn:]))
            e.ysens_post_exec = scipy.concatenate((e.ysens_post_exec, result[-1][eii].ysens_post_exec[n_dyn:]))
            e.ysens_pre_exec = scipy.concatenate((e.ysens_pre_exec, result[-1][eii].ysens_pre_exec[n_dyn:]))

def integrate_sensitivity(net, times, params=None, rtol=None, 
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

    # We use the master's result for events that occurred
    event_info = result[-2]
    events_occurred = result[-1]

    # Copy the sensitivity results into yout and (if necessary) youtdt
    # We don't need events_occurred here, because we already have it.
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
        if vars_assigned[worker]:
            _parse_sens_result(result, net, vars_assigned[worker], yout, youtdt, 
                               events_occurred)

    if exception_raised:
        raise exception_raised

    ddv_dpTrajectory = Trajectory_mod.Trajectory(net, is_sens=True, 
                                                 holds_dt=return_derivs)
    if return_derivs:
        yout = scipy.concatenate((yout, youtdt), axis=1)
    ddv_dpTrajectory.appendSensFromODEINT(tout, yout, holds_dt = return_derivs)
    ddv_dpTrajectory.events_occurred = events_occurred
    ddv_dpTrajectory.event_info = event_info

    net.trajectory = ddv_dpTrajectory

    return ddv_dpTrajectory

def dyn_var_fixed_point(net, dv0=None, with_logs=True, xtol=1e-6, time=0,
                        stability=False, fsolve_factor=100, 
                        maxfev=10000):
    """
    Return the dynamic variables values at the closest fixed point of the net.

    dv0  Initial guess for the fixed point. If not given, the current state
         of the net is used.
    with_logs   If True, the calculation is done in terms of logs of variables,
                so that they cannot be negative.
    xtol Tolerance to aim for.
    time Time to plug into equations.
    stability   If True, return the stability for the fixed point. -1 indicates
                stable node, +1 indicates unstable node, 0 indicates saddle
    fsolve_factor   'factor' argument for fsolve. For more information, see 
                    help(scipy.optimize.fsolve). Should be in range 0.1 to 100.
    maxfev      'maxfev' argument for fsolve. For more information, see 
                    help(scipy.optimize.fsolve). Should be an integer > 1.
    """
    net.compile()

    if dv0 is None:
        dv0 = scipy.array(net.getDynamicVarValues())
    else:
        dv0 = scipy.asarray(dv0)

    consts = net.constantVarValues

    zeros = scipy.zeros(len(dv0), scipy.float_)
    if with_logs:
        # We take the absolute value of dv0 to avoid problems from small 
        #  numerical noise negative values.
        if scipy.any(dv0 <= 0):
            logger.warning('Non-positive values in initial guess for fixed '
                           'point and with_logs = True. Rounding them up to '
                           'double_tiny. The most negative value was %g.' 
                           % min(dv0))
            dv0 = scipy.maximum(dv0, _double_tiny_)

        # XXX: Would like to replace these with C'd versions, if it's holding
        #      any of our users up.
        def func(logy):
            return net.res_function_logdv(time, logy, zeros, consts)
        def fprime(logy):
            y = scipy.exp(logy)
            y = scipy.maximum(y, _double_tiny_)
            return net.dres_dc_function(time, y, zeros, consts)

        x0 = scipy.log(dv0)
        # To transform sigma_x to sigma_log_x, we divide by x. We can set
        #  our sigma_log_x use to be the mean of what our xtol would yield
        #  for each chemical.
        xtol = scipy.mean(xtol/dv0)
    else:
        def func(y):
            return net.res_function(time, y, zeros, consts)
        def fprime(y):
            return net.dres_dc_function(time, y, zeros, consts)
        x0 = dv0

    try:
        dvFixed, infodict, ier, mesg =\
                scipy.optimize.fsolve(func, x0=x0.copy(), full_output=True,
                                      fprime=fprime,
                                      xtol=xtol, maxfev=maxfev,
                                      factor=fsolve_factor)
    except (scipy.optimize.minpack.error, ArithmeticError), X:
        raise FixedPointException(('Failure in fsolve.', X))

    tiny = _double_epsilon_

    if with_logs:
        dvFixed = scipy.exp(dvFixed)

    if ier != 1:
        if scipy.all(abs(dvFixed) < tiny) and not scipy.all(abs(x0) < 1e6*tiny):
            # This is the case where the answer is zero, and our initial guess
            # was reasonably large. In this case, the solver fails because
            # it's looking at a relative tolerance, but it's not really a
            # failure.
            pass
        else:
            raise FixedPointException(mesg, infodict)

    if not stability:
        return dvFixed
    else:
        jac = net.dres_dc_function(time, dvFixed, zeros, consts)
        u = scipy.linalg.eigvals(jac)
        u = scipy.real_if_close(u)
        if scipy.all(u < 0):
            stable = -1
        elif scipy.all(u > 0):
            stable = 1
        else:
            stable = 0
        return (dvFixed, stable)

def find_ypic_sens(y, yp, time, var_types, rtol, atol_for_sens, constants, 
                   net, opt_var, redirect_msgs=False):
    # On some systems, the f2py'd functions don't like len(constants)=0.
    if len(constants) == 0:
        constants = [0]
    var_types = scipy.asarray(var_types)
    y = scipy.asarray(y, scipy.float_)
    yp = scipy.asarray(yp, scipy.float_)
    atol_for_sens = scipy.asarray(atol_for_sens)

    N_dyn = len(var_types)
    y = copy.copy(y)
    yp = copy.copy(yp)

    # Find the initial conditions for the normal variables
    y_norm, yp_norm = find_ics(y[:N_dyn], yp[:N_dyn], time,
                               var_types, rtol, atol_for_sens[:N_dyn],
                               net.constantVarValues, net, 
                               redirect_msgs=redirect_msgs)
    # Copy the updated values into our y and yp arrays
    y[:N_dyn] = y_norm
    yp[:N_dyn] = yp_norm
    
    # Now we can solve for yp for all the *non-algebraic* sensitivity variables
    # Notice that they are simply the appropriate residual evaluated when
    #  yp for the sens variable is zero.
    yp[N_dyn:][var_types == 1] = 0
    res = net.sens_rhs(time, y, yp, constants)
    yp[N_dyn:][var_types == 1] = res[N_dyn:][var_types == 1]

    return yp

def find_ics(y, yp, time, var_types, rtol, atol, constants, net, 
             redirect_msgs=False):
    # We use this to find consistent sets of initial conditions for our
    #  integrations. (We don't let ddaskr do it, because it doesn't calculate
    #  values for d(alg_var)/dt, and we need them for sensitivity integration.)

    # On some systems, the f2py'd functions don't like len(constants)=0.
    if len(constants) == 0:
        constants = [0]
    var_types = scipy.asarray(var_types)
    atol = scipy.asarray(atol)
    rtol = scipy.asarray(rtol)
    # Note that we're copying y and yprime
    y = scipy.array(y, scipy.float_)
    yp = scipy.array(yp, scipy.float_)

    N_alg = scipy.sum(var_types == -1)

    dv_typ_vals = scipy.asarray([net.get_var_typical_val(id)
                                 for id in net.dynamicVars.keys()])

    if N_alg:
        # First we calculate a consistent set of algebraic variable values
        alg_vars_guess = y[var_types == -1]
        alg_typ_vals = dv_typ_vals[var_types == -1]
        possible_guesses = [alg_vars_guess, alg_typ_vals, 
                            scipy.ones(N_alg, scipy.float_)]
        redir = Utility.Redirector()
        if redirect_msgs:
            redir.start()
        try:
            for guess in possible_guesses:
                sln, infodict, ier, mesg = \
                        scipy.optimize.fsolve(net.alg_res_func, x0 = guess,
                                              xtol = min(rtol), 
                                              args = (y, time, constants),
                                              full_output=True)
                sln = scipy.atleast_1d(sln)
                final_residuals = net.alg_res_func(sln, y, time, constants)
                if not scipy.any(abs(final_residuals) 
                                 > abs(atol[var_types == -1])):
                    # This is success.
                    break
            else:
                message = ('Failed to calculate consistent algebraic values in '
                           'network %s.' % net.get_id())
                raise Utility.SloppyCellException(message)
        finally:
            messages = redir.stop()

        # Now plug those values into the current y
        y[var_types == -1] = sln

    # The non-algebraic variable yprimes come straight from the residuals
    yp_non_alg = net.res_function(time, y, y*0, constants)[var_types == 1]
    yp[var_types == 1] = yp_non_alg
                                              
    if not N_alg:
        return y, yp

    # Now we need to figure out yprime for the algebraic vars
    curr_alg_yp = yp[var_types == -1]
    ones_arr = scipy.ones(N_alg, scipy.float_)
    # We try a range of possible guesses. Note that this is really just
    # a linear system, so we continue to have difficulties with this part of
    # the calculation, or if it becomes a slow-down, we should consider doing
    # it by a linear solve, rather than using fsolve.
    possible_guesses = [curr_alg_yp, alg_typ_vals, 
                        ones_arr,
                        scipy.mean(abs(yp)) * ones_arr, 
                        -scipy.mean(abs(yp)) * ones_arr,
                        max(abs(yp)) * ones_arr, 
                        -max(abs(yp)) * ones_arr]
    if redirect_msgs:
        redir.start()
    try:
        for guess in possible_guesses:
            sln, infodict, ier, mesg = \
                    scipy.optimize.fsolve(net.alg_deriv_func, x0 = guess,
                                          xtol = min(rtol),
                                          args = (y, yp, time, constants),
                                          full_output=True)
            sln = scipy.atleast_1d(sln)
            final_residuals = net.alg_deriv_func(sln, y, yp, time, constants)
            if not scipy.any(abs(final_residuals) > abs(atol[var_types == -1])):
                break
        else:
            raise Utility.SloppyCellException('Failed to calculate alg var '\
                                              'derivatives in network %s.' 
                                              % net.get_id())
    finally:
        messages=redir.stop()
    sln = scipy.atleast_1d(sln)
    yp[var_types == -1] = sln

    return y, yp

def _reduce_times(yout, tout, times):
    jj = 0
    maxjj = len(tout)
    for ii, twanted in enumerate(times):
        while (tout[jj] != twanted) & (jj<maxjj-1):
            jj += 1
        yout[ii] = yout[jj]

    yout = yout[:len(times)]

    return yout
