"""
Methods for evaluating the dynamics of Network.
"""
__docformat__ = "restructuredtext en"

import copy
import sets

import scipy
import scipy.interpolate
import scipy.optimize

import logging
logger = logging.getLogger('ReactionNetworks.Dynamics')

import SloppyCell.lsodar as lsodar
odeintr = lsodar.odeintr
import SloppyCell.daskr as daskr
daeint = daskr.daeint
import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList
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


def integrate(net, times, params=None, rtol=1e-6, fill_traj=True,
              return_events=False, return_derivs=False,
              redirect_msgs=True):
    """
    Integrate a Network, returning a Trajectory.

    net            The Network to integrate.
    times          A sequence of times to include in the integration output.
    params         Parameters for the integration. If params=None, the current
                   parameters of net are used.
    rtol           Relative error tolerance for this integration.
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
    """
    times = scipy.array(times)
    rtol = min(rtol, global_rtol)
    net.compile()

    if params is not None:
        net.update_optimizable_vars(params)

    # If you ask for time = 0, we'll assume you want dynamic variable values
    #  reset.
    if times[0] == 0:
        net.resetDynamicVariables()

    if getattr(net, 'integrateWithLogs', False):
        def func(logDV, time):
            dv = scipy.exp(logDV)
            ddv_dt = net.get_ddv_dt(dv, time)
            return ddv_dt / dv
        Dfun = None
        atol = rtol
        rtol = 0
        def root_func(logDV, time):
            dv = scipy.exp(logDV)
            ret = net.root_func(dv, time)
            return ret
    else:
        func = net.get_ddv_dt
        Dfun = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))
        root_func = net.root_func
        ### This is messy but if we have a variable whose typical value is 1.0e6 say
        ### drops to almost zero during an integration, then we have to avoid
        ### integration errors making it significantly negative
        if rtol is not None:
            atol = [min(rtol*net.get_var_typical_val(id),1e-6) for id
                    in net.dynamicVars.keys()]

    IC = net.getDynamicVarValues()
    start = times[0]
    yout, tout = scipy.zeros((0, len(IC)), scipy.float_), []
    # time derviative of trajectory
    youtdt = scipy.zeros((0,len(IC)),scipy.float_)
    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events
    te, ye, ie = [], [], []
    pendingEvents = {}
    # After firing an event, we integrate for root_grace_t without looking for
    #  events to prevent finding the same one over and over again due to
    #  numerical imprecision
    root_grace_t = (times[-1] - times[0])/1e6
    event_just_fired = False

    exception_raised = None

    while start < times[-1]:
        if round(start, 12) in pendingEvents.keys():
            IC = net.executeEvent(pendingEvents[round(start, 12)], IC, start)
            del pendingEvents[round(start, 12)]

        # If an event just fired, we integrate for root_grace_t without looking
        #  for events, to prevent detecting the same event several times.
        if event_just_fired:
            if getattr(net, 'integrateWithLogs', False):
                IC = scipy.log(IC)
            next_requested = scipy.compress(times > start, times)[0]
            integrate_to = min(start + root_grace_t, next_requested)
            try:
                temp = odeintr(func, copy.copy(IC), [start, integrate_to],
                           Dfun = Dfun,
                           mxstep = 10000, rtol = rtol, atol = atol,
                           int_pts = fill_traj,
                           full_output = True, return_derivs = return_derivs,
                               redirect_msgs=redirect_msgs)
            except Utility.SloppyCellException, X:
                ### need to return as much of the trajectory as we have so far
                # since integration failed, return all the traj we got
                exception_raised = X
                print 'Exception in grace integration'
                break
            if temp[5] < 0:
                break

            # We don't append the last point, to prevent a needless 'event
            #  looking' duplication of times in the trajectory.
            # temp[5] is the derivative of the trajectory w.r.t. time
            if getattr(net, 'integrateWithLogs', False):
                yout = scipy.concatenate( (yout, scipy.exp(temp[0][:-1])) )
                if return_derivs :
                    youtdt = scipy.concatenate( (youtdt,scipy.exp(temp[0][:-1])*temp[5][:-1]) )
                start, IC = temp[1][-1], copy.copy(scipy.exp(temp[0][-1]))
            else:
                yout = scipy.concatenate((yout, temp[0][:-1]))
                if return_derivs :
                    youtdt = scipy.concatenate( (youtdt,temp[5][:-1]) )
                start, IC = temp[1][-1], copy.copy(temp[0][-1])
            tout.extend(temp[1][:-1])
            event_just_fired = False
            logger.debug('Finished event grace time integration.')

        # If we have pending events, only integrate until the next one.
        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())
            curTimes = scipy.compress((times > start) & (times < nextEventTime),
                                      times)
            curTimes = scipy.concatenate(([start], curTimes, [nextEventTime]))
        else:
            curTimes = scipy.compress(times > start, times)
            curTimes = scipy.concatenate(([start], curTimes))

        if getattr(net, 'integrateWithLogs', False):
            IC = scipy.log(IC)
        try:
            temp = odeintr(func, copy.copy(IC), curTimes,
                       root_func = root_func,
                       root_term = [True]*len(net.events),
                       Dfun = Dfun,
                       mxstep = 10000, rtol = rtol, atol = atol,
                       int_pts = fill_traj,
                       insert_events = True,
                       full_output = True, return_derivs = return_derivs,
                           redirect_msgs=redirect_msgs)
        except Utility.SloppyCellException, X:
            ### need to return as much of the trajectory as we have so far
            # since integration failed, return all the traj we got
            exception_raised = X
            global failed_args
            global failed_kwargs
            temp = X.args[1]
            failed_args = (func, IC, curTimes)
            failed_kwargs = {'root_func': root_func,
                             'root_term': [True]*len(net.events),
                             'Dfun': Dfun,
                             'mxstep': 10000, 
                             'rtol': rtol, 
                             'atol': atol,
                             'int_pts': fill_traj,
                             'insert_events': True,
                             'full_output': True, 
                             'return_derivs': return_derivs}

        if getattr(net, 'integrateWithLogs', False):
            yout = scipy.concatenate((yout, scipy.exp(temp[0])))
            if return_derivs :
                youtdt = scipy.concatenate( (youtdt, scipy.exp(temp[0])*temp[5] ) )
        else:
            yout = scipy.concatenate((yout, temp[0]))
            if return_derivs :
                youtdt = scipy.concatenate( (youtdt, temp[5] ) )

        tout.extend(temp[1])
        start, IC = tout[-1], copy.copy(yout[-1])

        if exception_raised:
            break

        # Check the directions of our root crossings to see whether events
        #  actually fired.
        for te_this, ye_this, ie_this in zip(temp[2], temp[3], temp[4]):
            if getattr(net, 'integrateWithLogs', False): 
                ye_this = scipy.exp(ye_this)
            root_derivs = net.root_func_dt(ye_this,
                                           net.get_ddv_dt(ye_this, te_this),
                                           te_this)
            if root_derivs[ie_this] > 0:
                te.append(te_this)
                ye.append(ye_this)
                ie.append(ie_this)
                event_just_fired = True

        # If an event fired
        if event_just_fired:
            event = net.events[ie[-1]]
            logger.debug('Event %s fired at time=%g in network %s.'
                         % (event.id, te[-1], net.id))
            delay = net.fireEvent(event, yout[-1], te[-1])
            pendingEvents[round(te[-1] + delay, 12)] = event

    if len(yout) and len(tout):
        net.updateVariablesFromDynamicVars(yout[-1], tout[-1])

    if not fill_traj and not exception_raised:
        yout = _reduce_times(yout, tout, times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, times)
        tout = times
    elif reduce_space :
        filtered_times = [tout[i] for i in range(0,len(tout),reduce_space)]
        filtered_times = scipy.sort( list( sets.Set(scipy.concatenate( (filtered_times,times) ) ) ) )
        # make sure we don't miss data times
        yout = _reduce_times(yout, tout, filtered_times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, filtered_times)
        tout = filtered_times

    if hasattr(net, 'kludge'):
        return tout, yout

    trajectory = Trajectory_mod.Trajectory(net,holds_dt=return_derivs)
    if return_derivs :
        yout = scipy.concatenate((yout,youtdt),axis=1) # join columnwise

    trajectory.appendFromODEINT(tout, yout, holds_dt=return_derivs)
    trajectory.add_event_info((te,ye,ie))
    net.trajectory = trajectory

    # raise exception if exited integration loop prematurely
    if exception_raised:
        logger.warn('Integration ended prematurely in network %s on node %i.'
                    % (net.id, my_rank))
        raise exception_raised

    return trajectory




def integrate_J(net, times, rtol, atol, params=None, fill_traj=True,
              return_events=False, return_derivs=False,
              redirect_msgs=True):
    """
    Integrate a Network, returning a Trajectory.

    net            The Network to integrate.
    times          A sequence of times to include in the integration output.
    params         Parameters for the integration. If params=None, the current
                   parameters of net are used.
    rtol           Relative error tolerance for this integration.
    rtol           Absolute error tolerance for this integration.
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
    """

    # define a residual function that will work with the f2py interface
    def res_func(time, y, yprime, ires):
        return net.res_function(time, y, yprime, ires) 

    times = scipy.array(times)
    net.compile()

    if params is not None:
        net.update_optimizable_vars(params)

    # If you ask for time = 0, we'll assume you want dynamic variable values
    #  reset.
    if times[0] == 0:
        net.resetDynamicVariables()


    # removed 'integrateWithLogs' functionality for now
    """

    if getattr(net, 'integrateWithLogs', False):
        def func(logDV, time):
            dv = scipy.exp(logDV)
            ddv_dt = net.get_ddv_dt(dv, time)
            return ddv_dt / dv
        Dfun = None
        atol = rtol
        rtol = 0
        def root_func(logDV, time):
            dv = scipy.exp(logDV)
            ret = net.root_func(dv, time)
            return ret
    else:
        func = net.get_ddv_dt
        Dfun = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))
        root_func = net.root_func
        ### This is messy but if we have a variable whose typical value is 1.0e6 say
        ### drops to almost zero during an integration, then we have to avoid
        ### integration errors making it significantly negative
        if rtol is not None:
            atol = [min(rtol*net.get_var_typical_val(id),1e-6) for id
                    in net.dynamicVars.keys()]
    """

    # get this initial state, time, and derivative.
    IC = net.getDynamicVarValues()
    start = times[0]
    youtdt = net.get_ddv_dt(IC,start)

    # start variables for output storage
    yout, tout = scipy.zeros((0, len(IC)), scipy.float_), []

    # Will eventually need to fix this Dfun
    # Dfun = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))

    root_func = net.ddaskr_root
    num_events = len(net.events)

    
    ### This is messy but if we have a variable whose typical value is 1.0e6 say
    ### drops to almost zero during an integration, then we have to avoid
    ### integration errors making it significantly negative
    # out for now
    """
    if rtol is not None:
        atol = [min(rtol*net.get_var_typical_val(id),1e-6) for id
                in net.dynamicVars.keys()]
    """
    
    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events
    te, ye, ie = [], [], []
    pendingEvents = {}

    # After firing an event, we integrate for root_grace_t without looking for
    #  events to prevent finding the same one over and over again due to
    #  numerical imprecision
    root_grace_t = (times[-1] - times[0])/1e6
    event_just_fired = False

    exception_raised = None
    
    while start < times[-1]:

        # not doing delays right now
        if round(start, 12) in pendingEvents.keys():
            IC = net.executeEvent(pendingEvents[round(start, 12)], IC, start)
            del pendingEvents[round(start, 12)]


        # If an event just fired, we integrate for root_grace_t without looking
        #  for events, to prevent detecting the same event several times.
        if event_just_fired:
            if getattr(net, 'integrateWithLogs', False):
                IC = scipy.log(IC)
            next_requested = scipy.compress(times > start, times)[0]
            integrate_to = min(start + root_grace_t, next_requested)
            try:
                temp = daeint(res_func, [start, integrate_to],
                              copy.copy(IC), copy.copy(youtdt),
                              rtol, atol,
                              intermediate_output = fill_traj,
                              max_steps = 10000,
                              redir_output = redirect_msgs)

            except Utility.SloppyCellException, X:
                ### need to return as much of the trajectory as we have so far
                # since integration failed, return all the traj we got
                exception_raised = X
                print 'Exception in grace integration'
                break
            if temp[5] < 0:
                break

            # We don't append the last point, to prevent a needless 'event
            #  looking' duplication of times in the trajectory.
            # temp[5] is the derivative of the trajectory w.r.t. time
            if getattr(net, 'integrateWithLogs', False):
                yout = scipy.concatenate( (yout, scipy.exp(temp[0][:-1])) )
                if return_derivs :
                    youtdt = scipy.concatenate( (youtdt,scipy.exp(temp[0][:-1])*temp[5][:-1]) )
                start, IC = temp[1][-1], copy.copy(scipy.exp(temp[0][-1]))
            else:
                yout = scipy.concatenate((yout, temp[0][:-1]))
                if return_derivs :
                    youtdt = scipy.concatenate( (youtdt,temp[5][:-1]) )
                start, IC = temp[1][-1], copy.copy(temp[0][-1])
            tout.extend(temp[1][:-1])
            event_just_fired = False
            logger.debug('Finished event grace time integration.')


        # If we have pending events, only integrate until the next one.
        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())
            curTimes = scipy.compress((times > start) & (times < nextEventTime),
                                      times)
            curTimes = scipy.concatenate(([start], curTimes, [nextEventTime]))
        else:
            curTimes = scipy.compress(times > start, times)
            curTimes = scipy.concatenate(([start], curTimes))
        """
        if getattr(net, 'integrateWithLogs', False):
            IC = scipy.log(IC)
        """

        try:
            temp = daeint(res_func, curTimes,
                          copy.copy(IC), copy.copy(youtdt),
                          rtol, atol,
                          rt = root_func,
                          intermediate_output = fill_traj,
                          redir_output = redirect_msgs)            
            
        except Utility.SloppyCellException, X:
            ### need to return as much of the trajectory as we have so far
            # since integration failed, return all the traj we got
            exception_raised = X
            global failed_args
            global failed_kwargs
            temp = X.args[1]
            failed_args = (func, IC, curTimes)
            failed_kwargs = {'root_func': root_func,
                             'root_term': [True]*len(net.events),
                             'Dfun': Dfun,
                             'mxstep': 10000, 
                             'rtol': rtol, 
                             'atol': atol,
                             'int_pts': fill_traj,
                             'insert_events': True,
                             'full_output': True, 
                             'return_derivs': return_derivs}



        if getattr(net, 'integrateWithLogs', False):
            yout = scipy.concatenate((yout, scipy.exp(temp[0])))
            if return_derivs :
                youtdt = scipy.concatenate( (youtdt, scipy.exp(temp[0])*temp[5] ) )
        else:
            yout = scipy.concatenate((yout, temp[0]))
            if return_derivs :
                youtdt = scipy.concatenate( (youtdt, temp[2] ) )
            

        tout.extend(temp[1])
        start, IC = tout[-1], copy.copy(yout[-1])

        if exception_raised:
            break

        # Check the directions of our root crossings to see whether events
        # actually fired.  DASKR automatically returns the direction of event
        # crossings, so we can read this information from the integration results
        # zip(...) allows us to iterate over multiple sequences in parallel.
        for ii, dire_this in zip(range(num_events), temp[5]):
            """
            if getattr(net, 'integrateWithLogs', False): 
                ye_this = scipy.exp(ye_this)
            """
            
            # dire_this is the direction of the event crossing.  If it's
            # 1 that means Ri changed from negative to positive.  In that case
            # we need to record the event.
        
            if dire_this > 0:
                # an event has fired, so record the state and event number
                te.append(temp[3])
                ye.append(temp[4])
                ie.append(ii)
                event_just_fired = True

        # If an event fired
        if event_just_fired:
            event = net.events[ie[-1]]
            logger.debug('Event %s fired at time=%g in network %s.'
                         % (event.id, te[-1], net.id))
            delay = net.fireEvent(event, yout[-1], te[-1])
            pendingEvents[round(te[-1] + delay, 12)] = event


    if len(yout) and len(tout):
        net.updateVariablesFromDynamicVars(yout[-1], tout[-1])

    if not fill_traj and not exception_raised:
        yout = _reduce_times(yout, tout, times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, times)
        tout = times
    elif reduce_space :
        filtered_times = [tout[i] for i in range(0,len(tout),reduce_space)]
        filtered_times = scipy.sort( list( sets.Set(scipy.concatenate( (filtered_times,times) ) ) ) )
        # make sure we don't miss data times
        yout = _reduce_times(yout, tout, filtered_times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, filtered_times)
        tout = filtered_times

    if hasattr(net, 'kludge'):
        return tout, yout

    trajectory = Trajectory_mod.Trajectory(net,holds_dt=return_derivs)
    if return_derivs :
        yout = scipy.concatenate((yout,youtdt),axis=1) # join columnwise

    trajectory.appendFromODEINT(tout, yout, holds_dt=return_derivs)
    trajectory.add_event_info((te,ye,ie))
    net.trajectory = trajectory

    # raise exception if exited integration loop prematurely
    if exception_raised:
        logger.warn('Integration ended prematurely in network %s on node %i.'
                    % (net.id, my_rank))
        raise exception_raised

    return trajectory


def integrate_sens_subset(net, times, rtol=1e-6, 
                          fill_traj=False, opt_vars=None, 
                          return_derivs=False):
    nDyn, nOpt = len(net.dynamicVars), len(opt_vars)

    atolDv = [min(rtol*net.get_var_typical_val(id), 1e-6) for id
              in net.dynamicVars.keys()]

    start = times[0]
    # IC is zero unless IC is a function of parameters
    IC = scipy.zeros(nDyn * (nOpt + 1), scipy.float_)

    # Handle sensitivities of initial conditions
    for dvInd, (id, var) in enumerate(net.dynamicVars.items()):
        if isinstance(var.initialValue, str):
            for ovInd, ovName in enumerate(opt_vars) :
                DwrtOV = net.takeDerivative(var.initialValue, ovName)
                IC[nDyn*(ovInd+1) + dvInd] = net.evaluate_expr(DwrtOV,
                                                            time=start)

    # Copy in basic ICs
    IC[:nDyn] = net.getDynamicVarValues()

    yout, tout = scipy.zeros((0, len(IC)), scipy.float_), []
    youtdt = scipy.zeros((0, len(IC)), scipy.float_)

    te, ye, ie = [], [], []
    pendingEvents = {}

    while start < times[-1]:
        if start in pendingEvents.keys():
            IC = net.executeEvent(pendingEvents[start], IC, start)
            pendingEvents.pop(start)

        if pendingEvents:
            nextEventTime = min(pendingEvents.keys())
            curTimes = scipy.compress((times > start) & (times < nextEventTime),
                                    times)
            curTimes = scipy.concatenate(([start], curTimes, [nextEventTime]))
        else:
            curTimes = scipy.compress(times > start, times)
            curTimes = scipy.concatenate(([start], curTimes))

        # We'll do an initial integration run without the sensitivites
        #  to find out how far we need to integrate to get to a terminating
        #  event. This reduces the amount of backtracking we need to do
        #  when doing the much slower sensitivity integrations.
        func = net.get_ddv_dt
        nRoots = len(net.events)
        Dfun_temp = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))

        try:
            temp = odeintr(func, copy.copy(IC[:nDyn]),
                           curTimes,
                           root_func = net.root_func,
                           root_term = [True]*nRoots,
                           Dfun = Dfun_temp,
                           mxstep = 10000, rtol = rtol, atol = atolDv,
                           int_pts = fill_traj,
                           insert_events = True,
                           full_output = True, return_derivs = return_derivs)
        except Utility.SloppyCellException, X:
            logger.warn('Sensitivity integration failed for network %s on '
                        'node %i.' % (net.id, my_rank))
            raise

        curTimes = temp[1]
        newSpace = scipy.zeros((len(curTimes), nDyn * (nOpt + 1)), scipy.float_)
        yout = scipy.concatenate((yout, newSpace))
        yout[-len(curTimes):, :nDyn] = copy.copy(temp[0][:,:nDyn])
        if return_derivs :
            youtdt = scipy.concatenate((youtdt, newSpace))
            youtdt[-len(curTimes):, :nDyn] = copy.copy(temp[5][:,:nDyn])

        tout.extend(curTimes)
        te.extend(temp[2])
        ye.extend(temp[3])
        ie.extend(temp[4])

        for ovIndex, ovId in enumerate(opt_vars):
            ICTrunc = scipy.zeros(2 * nDyn, scipy.float_)
            ICTrunc[:nDyn] = copy.copy(IC[:nDyn])
            ICTrunc[nDyn:] = copy.copy(IC[(ovIndex + 1) * nDyn:
                                        (ovIndex + 2) * nDyn])

            # Absolute tolerances for sensitivity variables
            if rtol is None or scipy.isscalar(rtol):
                rtolForThis = rtol
            else:
                rtolForThis = scipy.concatenate((rtol, rtol))
            if atolDv is not None:
                atolSens = scipy.asarray(atolDv)/net.get_var_typical_val(ovId)
                atolForThis = scipy.concatenate((atolDv, atolSens))
            else:
                atolForThis = None

            try:
                temp2 = odeintr(_Ddv_and_DdvDov_dtTrunc,
                                ICTrunc, curTimes,
                                args = (net, 
                                        net.optimizableVars.indexByKey(ovId)),
                                mxstep = 10000,
                                rtol = rtolForThis, atol = atolForThis, 
                                return_derivs = return_derivs)
            except Utility.SloppyCellException, X:
                logger.warn('Sensitivity integration failed for network %s on '
                            'node %i during optimizable variable %s.'
                            % (net.id, my_rank, ovId))
                raise


            yout[-len(curTimes):, (ovIndex + 1) * nDyn:(ovIndex + 2)*nDyn] =\
                    copy.copy(temp2[0][:,nDyn:])

            # temp2[2] is the derivative of the trajectory w.r.t time now 
            # because no events are being returned
            if return_derivs :
                youtdt[-len(curTimes):, (ovIndex + 1) * nDyn:(ovIndex + 2)*nDyn] =\
                        copy.copy(temp2[5][:,nDyn:])

        start, IC = tout[-1], copy.copy(yout[-1])

        if temp[4]:
            event = net.events[ie[-1]]
            delay = net.fireEventAndUpdateDdv_Dov(event, yout[-1], te[-1], opt_vars)
            pendingEvents[te[-1] + delay] = event

    # If we broke_out, fill the trajectory with what we have
    if not fill_traj and len(yout) and len(tout):
        yout = _reduce_times(yout, tout, times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, times)
        tout = times
    elif reduce_space :
        filtered_times = [tout[i] for i in range(0,len(tout),reduce_space)]
        filtered_times = scipy.sort( list( sets.Set(scipy.concatenate( (filtered_times,times) ) ) ) )
        # make sure we don't miss data times
        yout = _reduce_times(yout, tout, filtered_times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, filtered_times)
        tout = filtered_times

    if not return_derivs:
        return tout, yout
    elif return_derivs:
        return tout, yout, youtdt

def _parse_sens_result(result, net, opt_vars, yout, youtdt=None):
    """
    Utility function for parsing the return from integrate_sens_subset
    """
    n_dyn, n_opt = len(net.dynamicVars), len(net.optimizableVars)
    for work_index, opt_var in enumerate(opt_vars):
        net_index = net.optimizableVars.indexByKey(opt_var)
        yout[:, (net_index+1)*n_dyn: (net_index+2)*n_dyn] =\
                result[1][:, (work_index+1)*n_dyn:(work_index+2)*n_dyn]
        if youtdt is not None:
            youtdt[:, (net_index+1)*n_dyn: (net_index+2)*n_dyn] =\
                    result[2][:, (work_index+1)*n_dyn:(work_index+2)*n_dyn]


def integrate_sensitivity(net, times, params=None, rtol=1e-6, 
                          fill_traj=False, return_derivs=False):
    logger.debug('Entering integrate_sens on node %i' % my_rank)
    if not scipy.isscalar(rtol):
        rtol = 1e-6

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
                'rtol, fill_traj, opt_vars, return_derivs)'
        args = {'net':net, 'times':times, 'rtol':rtol, 'fill_traj':fill_traj,
                'opt_vars':vars_assigned[worker], 'return_derivs':return_derivs}
        pypar.send((command, args), worker)

    logger.debug('Master doing vars %s' % str(vars_assigned[0]))
    try:
        result = integrate_sens_subset(net, times, rtol, fill_traj, 
                                       vars_assigned[0], return_derivs)
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
        dvFixed, infodict, ier, mesg = scipy.optimize.fsolve(func, x0=x0,
                                                             fprime=fprime,
                                                             col_deriv=True,
                                                             full_output=True,
                                                             xtol=xtol)
    except scipy.optimize.minpack.error, X:
        raise FixedPointException(('Failure in fsolve.', X))

    if ier != 1:
        raise FixedPointException(mesg, infodict)

    if with_logs:
        return scipy.exp(dvFixed)
    else:
        return dvFixed

def _Ddv_and_DdvDov_dtTrunc_2(sens_y, time, net, ovIndex, tcks):
    dv_y = [scipy.interpolate.splev(time, tck) for tck in tcks]
    D2dv_Dov_dt = _D2dv_Dov_dtTrunc(dv_y, sens_y, time, net, ovIndex)

    return D2dv_Dov_dt

def Dfun_sens(sens_y, time, net, ovIndex, tcks):
    dv_y = [scipy.interpolate.splev(time, tck) for tck in tcks]
    return scipy.transpose(net.get_d2dv_ddvdt(dv_y,time))

def _D2dv_Dov_dtTrunc(dv, Ddv_Dov, time, net, ovIndex):
    # The partial derivative of the dynamic vars wrt the ovIndex optimizable var
    partials = net.get_d2dv_dovdt(dv, time, indices = [ovIndex])[:, ovIndex]

    # The partial derivative of the dynamic vars wrt the other dynamic vars
    J = net.get_d2dv_ddvdt(dv, time)

    # Now we do the chain rule to the the total derivative of the dynamic vars
    #  wrt the specified optimizable var.
    # Transpose of J is [[dAdA, dAdB],[dBdA, dBdB]] and Ddv_dov is [dAda, dBda]
    #  This product is [dAdA*dAda + dAdB*dBdb, dBdA*dAda + dBdB*dBdb]
    #  (in a notation where A is a dynamic var and a is an optimizable)
    D2dv_Dov_dt = scipy.dot(scipy.transpose(J), Ddv_Dov) + partials

    return D2dv_Dov_dt

def _Ddv_and_DdvDov_dtTrunc(y, time, net, ovIndex):
    nDV = len(net.dynamicVars)

    # XXX: Might want to fill in a global array here. Perhaps a small speed
    #      boost
    Dc_Dt = net.get_ddv_dt(y[:nDV], time)
    D2dv_Dov_dt = _D2dv_Dov_dtTrunc(y[:nDV], y[nDV:], time, net, ovIndex)

    return scipy.concatenate((Dc_Dt, D2dv_Dov_dt))

def _reduce_times(yout, tout, times):
    jj = 0
    maxjj = len(tout)
    for ii, twanted in enumerate(times):
        while (tout[jj] != twanted) & (jj<maxjj-1):
            jj += 1
        yout[ii] = yout[jj]

    yout = yout[:len(times)]

    return yout


def integrate_sensitivity_2(net, times, params=None, rtol = 1e-6):
    if HAVE_PYPAR:
        import pypar
        numproc = pypar.size()
        myid = pypar.rank()
        node = pypar.get_processor_name()
    else:
        numproc, myid, node = 1, 0, 'serial'

    rtol = min(rtol, global_rtol)
    traj = integrate(net, times, params, fill_traj=True,
                     return_events=True, rtol=rtol/10.)

    te, ye, ie = traj.event_info # event info stored in trajectory

    # should we distinguish between passed in times and these times?
    # otherwise we always have a filled in trajectory
    times = traj.get_times()
    n_dyn, n_opt = len(net.dynamicVars), len(net.optimizableVars)

    start = times[0]
    # IC is zero unless IC is a function of parameters
    IC = scipy.zeros(n_dyn * (n_opt + 1), scipy.float_)
    # In the case where return_derivs = True, start_vals still only contains dynamic
    # vars and assigned vars, because the time derivatives are not part of the network's
    # variable list (which Trajectory_mod uses)
    start_vals = traj.get_var_vals(start)
    net.set_var_vals(start_vals)
    # Handle sensitivities of initial conditions
    for dvInd, id in enumerate(net.dynamicVars.keys()):
        init_val = net.get_var_ic(id)
    if isinstance(init_val, str):
            for ovInd, ovName in enumerate(net.optimizableVars.keys()) :
                DwrtOV = net.takeDerivative(init_val, ovName)
                IC[n_dyn*(ovInd+1) + dvInd] = net.evaluate_expr(DwrtOV,
                                                               time=start)

    yout = scipy.zeros((len(times), len(IC)), scipy.float_)
    youtdt = scipy.zeros((len(times), len(IC)), scipy.float_)

    for ii, dyn_id in enumerate(net.dynamicVars.keys()):
        yout[:,ii] = traj.get_var_traj(dyn_id)
        if return_derivs:
            youtdt[:,ii] = traj.get_var_traj((dyn_id,'time'))
    pendingEvents = {}

    intervals = computeIntervals(traj)
    dv_typ_vals = [net.get_var_typical_val(dyn_id) for dyn_id
                   in net.dynamicVars.keys()]

    # Dish out work to processes if running in paralell
    load = [int(scipy.floor(n_opt/numproc))]*numproc
    # now deal with the remainder -- spread it out
    for i in range(n_opt%numproc) :
        load[i] = load[i] + 1
    cur_ovIndex = [0]*numproc
    last_ovIndex = [0]*numproc

    for curid in range(numproc) :
        for i in range(curid) :
            cur_ovIndex[curid] += load[i]

    all_ovs = [ (ovIndex, ovId) for ovIndex,ovId in enumerate(net.optimizableVars.keys())]

    for curid in range(numproc) :
        last_ovIndex[curid] = cur_ovIndex[curid] + load[curid] - 1


    for start_ind, end_ind in intervals:
        start_time = times[start_ind]
        end_time = times[end_ind-1]

        if start_time in pendingEvents.keys():
            # We don't currently do delays in sensitivity integration, but
            #  we'll keep this here for when we do.
            IC = net.executeEvent(pendingEvents[start_time], IC, start_time)
            pendingEvents.pop(start_time)

        curTimes = traj.get_times()[start_ind:end_ind]

        # The order of the spline to fit. If we have too few points, we have
        #  to go smaller than 5.
        k = min(5, end_ind - start_ind - 1)
        ys = [traj.get_var_traj(dv_id)[start_ind:end_ind] for dv_id
              in net.dynamicVars.keys()]
        tcks = [scipy.interpolate.splrep(curTimes, y, k=k, s=0) for y in ys]

        for ov_ind, opt_id in all_ovs[cur_ovIndex[myid]:(last_ovIndex[myid]+1)] :
            IC_this= copy.copy(IC[(ov_ind + 1) * n_dyn:
                                   (ov_ind + 2) * n_dyn])

            atol = rtol * scipy.asarray(dv_typ_vals)/\
                    net.get_var_typical_val(opt_id)


            temp = odeintr(_Ddv_and_DdvDov_dtTrunc_2,
                           IC_this, curTimes,
                           args = (net, ov_ind, tcks),
                           Dfun = Dfun_sens,
                           mxstep = 10000,
                           rtol = rtol, atol = atol, return_derivs = return_derivs,
                           tcrit = (times[-1],))

            yout[start_ind:end_ind, (ov_ind + 1) * n_dyn:(ov_ind + 2)*n_dyn] =\
                    copy.copy(temp[0])
            if return_derivs :
                youtdt[start_ind:end_ind,(ov_ind + 1) * n_dyn:(ov_ind + 2)*n_dyn] =\
                    copy.copy(temp[5])

        IC = copy.copy(yout[end_ind - 1])

        if end_time in te:
            ind = te.index(end_time)
            event = net.events[ie[ind]]
            delay = net.fireEventAndUpdateDdv_Dov(event, yout[-1], te[ind])
            pendingEvents[te[ind] + delay] = event

    net.updateVariablesFromDynamicVars(yout[-1][:n_dyn], times[-1])

    tout = times[:] # copy times
    if not net.add_int_times: # this doesn't mean anything because times passed in is overwritten
        yout = _reduce_times(yout, tout, times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, times)
        tout = times
    elif reduce_space :
        filtered_times = [tout[i] for i in range(0,len(tout),reduce_space)]
        filtered_times = scipy.sort( list( sets.Set(scipy.concatenate( (filtered_times,times) ) ) ) )
        yout = _reduce_times(yout, tout, filtered_times)
        if return_derivs :
            youtdt = _reduce_times(youtdt, tout, filtered_times)
        tout = filtered_times

    ddv_dpTrajectory = Trajectory_mod.Trajectory(net, is_sens=True, holds_dt = return_derivs)

    if myid==0 and numproc>1: # I am the master
        for procno in range(1,numproc) :
            logger.debug('%s receiving package from %d ' % (procno,pypar.Get_processor_name()) )
            yout[:,(cur_ovIndex[procno]+1)*n_dyn:(last_ovIndex[procno]+2)*n_dyn] = pypar.receive(procno,tag=1)
            if return_derivs :
                youtdt[:,(cur_ovIndex[procno]+1)*n_dyn:(last_ovIndex[procno]+2)*n_dyn] = pypar.receive(procno,tag=2)

    elif numproc>1 :  # I am a slave
        logger.debug('%s sending package ' % pypar.Get_processor_name())
        pypar.send(yout[:,(cur_ovIndex[myid] + 1)*n_dyn:(last_ovIndex[myid] + 2)*n_dyn],0,tag=1)  # send to the master
        if return_derivs :
            pypar.send(youtdt[:,(cur_ovIndex[myid] + 1)*n_dyn:(last_ovIndex[myid] + 2)*n_dyn],0,tag=2)

    if return_derivs :
        yout = scipy.concatenate((yout,youtdt),axis=1)
    ddv_dpTrajectory.appendSensFromODEINT(tout, yout, holds_dt=return_derivs)

    return ddv_dpTrajectory

def computeIntervals(traj):
    # We want to break up our integrals when events fire, so first we figure out
    #  when they fired by looking for duplicated times in the trajectory
    eventIndices = scipy.compress(scipy.diff(traj.get_times()) == 0,
                                  scipy.arange(len(traj.get_times())))
    intervals = zip([0] + list(eventIndices + 1),
                    list(eventIndices + 1) + [len(traj.timepoints)])

    return intervals

try:
    import psyco
    psyco.bind(scipy.interpolate.splrep)
    psyco.bind(scipy.interpolate.splev)
    psyco.bind(_D2dv_Dov_dtTrunc)
    psyco.bind(_Ddv_and_DdvDov_dtTrunc)
    psyco.bind(_Ddv_and_DdvDov_dtTrunc_2)
    psyco.bind(SloppyCell.lsodar.odeintr)
    psyco.bind(Dfun_sens)
except ImportError:
    pass
