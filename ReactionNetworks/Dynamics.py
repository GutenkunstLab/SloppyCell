import copy
import sets

import scipy

import SloppyCell.lsodar
odeintr = SloppyCell.lsodar.odeintr
import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList
import Trajectory_mod

# XXX: Need to incorporate root_grace_t into integrate_sensitivity
global_rtol = 1e-6

def integrate(net, times, params=None, rtol=1e-6, fill_traj=False,
              return_events=False):
    rtol = min(rtol, global_rtol)
    net.compile()
    fill_traj = net.add_int_times or fill_traj

    if params is not None:
        net.update_optimizable_vars(params)

    # Do we want to add times to the tail of our integration for pretty plotting
    #  with data?
    times = scipy.array(times)
    if net.add_tail_times:
        times = scipy.concatenate((times, [1.05*times[-1]]))

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
    else:
        func = net.get_ddv_dt
        Dfun = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))
        if rtol is not None:
            atol = [rtol*net.get_var_typical_val(id) for id
                    in net.dynamicVars.keys()]

    IC = net.getDynamicVarValues()
    start = times[0]
    yout, tout = scipy.zeros((0, len(IC)), scipy.Float), []
    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events
    te, ye, ie = [], [], []
    pendingEvents = {}
    root_grace_t = (times[-1] - times[0])/1e10
    event_just_fired = False

    while start < times[-1]:
        if round(start, 12) in pendingEvents.keys():
            IC = net.executeEvent(pendingEvents[round(start, 12)], IC, start)
            del pendingEvents[round(start, 12)]

        # If an event just fired, we integrate for root_grace_t without looking
        #  for events, to prevent detecting the same event several times.
        if event_just_fired:
            if getattr(net, 'integrateWithLogs', False):
                IC = scipy.log(IC)
            temp = odeintr(func, copy.copy(IC), [start, start+root_grace_t], 
                           Dfun = Dfun, 
                           mxstep = 10000, rtol = rtol, atol = atol,
                           int_pts = fill_traj,
                           full_output = True)

            if getattr(net, 'integrateWithLogs', False):
                yout = scipy.concatenate((yout, scipy.exp(temp[0])))
            else:
                yout = scipy.concatenate((yout, temp[0]))
            tout.extend(temp[1])
            start, IC = tout[-1], copy.copy(yout[-1])
            event_just_fired = False

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
        temp = odeintr(func, copy.copy(IC), curTimes, 
                       root_func = net.root_func, 
                       root_term = [True]*len(net.events),
                       Dfun = Dfun, 
                       mxstep = 10000, rtol = rtol, atol = atol,
                       int_pts = fill_traj,
                       insert_events = True,
                       full_output = True)

        if getattr(net, 'integrateWithLogs', False):
            yout = scipy.concatenate((yout, scipy.exp(temp[0])))
        else:
            yout = scipy.concatenate((yout, temp[0]))
        tout.extend(temp[1])
        start, IC = tout[-1], copy.copy(yout[-1])

        # Check the directions of our root crossings to see whether events
        #  actually fired.
        for te_this, ye_this, ie_this in zip(temp[2], temp[3], temp[4]):
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
            delay = net.fireEvent(event, yout[-1], te[-1])
            pendingEvents[round(te[-1] + delay, 12)] = event

    net.updateVariablesFromDynamicVars(yout[-1], tout[-1])

    if not fill_traj:
        yout = _reduce_times(yout, tout, times)
        tout = times

    trajectory = Trajectory_mod.Trajectory(net)
    trajectory.appendFromODEINT(tout, yout)

    if return_events:
        return trajectory, te, ye, ie
    else:
        return trajectory

def integrate_sensitivity(net, times, params=None, rtol=1e-6):
    net.compile()

    if params is not None:
        net.update_optimizable_vars(params)

    times = scipy.array(times)
    if net.add_tail_times:
        times = scipy.concatenate((times, [1.05*times[-1]]))

    # If you ask for time = 0, we'll assume you want dynamic variable values
    #  reset.
    if times[0] == 0:
        net.resetDynamicVariables()

    nDyn, nOpt = len(net.dynamicVars), len(net.optimizableVars)

    atolDv = None
    if scipy.isscalar(rtol):
        atolDv = [rtol*net.get_var_typical_val(id) for id
                  in net.dynamicVars.keys()]
    elif rtol is not None:
        if len(rtol) == len(IC):
            atol = [rtol[ii]*net.get_var_typical_val(id) for (ii, id)
                    in enumerate(net.dynamicVars.keys())]
        else:
            sys.stderr.write("Non-scalar rtol of improper length passed into integrate_sensitivity!")

    start = times[0]
    # IC is zero unless IC is a function of parameters
    IC = scipy.zeros(nDyn * (nOpt + 1), scipy.Float) 

    # Handle sensitivities of initial conditions
    for dvInd, (id, var) in enumerate(net.dynamicVars.items()):
	if isinstance(var.initialValue, str):
            for ovInd, ovName in enumerate(net.optimizableVars.keys()) :
                DwrtOV = net.takeDerivative(var.initialValue, ovName)
                IC[nDyn*(ovInd+1) + dvInd] = net.evaluate_expr(DwrtOV,
                                                               time=start)

    # Copy in basic ICs
    IC[:nDyn] = net.getDynamicVarValues()

    #yout, tout = scipy.array([IC]), [start]
    yout, tout = scipy.zeros((0, len(IC)), scipy.Float), []
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

        temp = odeintr(func, copy.copy(IC[:nDyn]),
                       curTimes,
                       root_func = net.root_func,
                       root_term = [True]*nRoots,
                       Dfun = Dfun_temp,
                       mxstep = 10000, rtol = rtol, atol = atolDv,
                       int_pts = net.add_int_times,
                       insert_events = True,
                       full_output = True)

        curTimes = temp[1]
        newSpace = scipy.zeros((len(curTimes), nDyn * (nOpt + 1)), scipy.Float)
        yout = scipy.concatenate((yout, newSpace))
        yout[-len(curTimes):, :nDyn] = copy.copy(temp[0][:,:nDyn])
        tout.extend(curTimes)
        te.extend(temp[2])
        ye.extend(temp[3])
        ie.extend(temp[4])

        for ovIndex, ovId in enumerate(net.optimizableVars.keys()):
            ICTrunc = scipy.zeros(2 * nDyn, scipy.Float)
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

            temp2 = odeintr(_Ddv_and_DdvDov_dtTrunc,
                            ICTrunc, curTimes,
                            args = (net, ovIndex),
                            mxstep = 10000, 
                            rtol = rtolForThis, atol = atolForThis)

            yout[-len(curTimes):, (ovIndex + 1) * nDyn:(ovIndex + 2)*nDyn] =\
                    copy.copy(temp2[0][:,nDyn:])

        start, IC = tout[-1], copy.copy(yout[-1])

        if temp[4]:
            event = net.events[ie[-1]]
            delay = net.fireEventAndUpdateDdv_Dov(event, yout[-1], te[-1])
            pendingEvents[te[-1] + delay] = event

    net.updateVariablesFromDynamicVars(yout[-1][:nDyn], tout[-1])

    if not net.add_int_times:
        yout = _reduce_times(yout, tout, times)
        tout = times

    ddv_dpTrajectory = Trajectory_mod.Trajectory(net, is_sens=True)
    ddv_dpTrajectory.appendSensFromODEINT(tout, yout)

    return ddv_dpTrajectory

def dyn_var_fixed_point(net, dv0 = None):
    if dv0 is None:
        dv0 = net.getDynamicVarValues()

    ddv_dtFromLogs = lambda logDV: net.get_ddv_dt(scipy.exp(logDV), 0)
    fprime = lambda logDV: net.get_d2dv_ddvdt(scipy.exp(logDV), 0)\
            *scipy.exp(logDV)

    dvFixed = scipy.optimize.fsolve(ddv_dtFromLogs, x0 = scipy.log(dv0),
                                    fprime = fprime, col_deriv = True)

    return scipy.exp(dvFixed)

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
    for ii, twanted in enumerate(times):
        while tout[jj] != twanted:
            jj += 1
        yout[ii] = yout[jj]

    yout = yout[:len(times)]

    return yout

def integrate_sensitivity_2(net, times, params=None, rtol = 1e-6):
    rtol = min(rtol, global_rtol)
    traj, te, ye, ie = integrate(net, times, params, fill_traj=True,
                                 return_events=True, rtol=1e-7)

    times = traj.get_times()
    n_dyn, n_opt = len(net.dynamicVars), len(net.optimizableVars)
    IC = scipy.zeros(n_dyn * (n_opt + 1), scipy.Float)

    start = times[0]
    # IC is zero unless IC is a function of parameters
    IC = scipy.zeros(n_dyn * (n_opt + 1), scipy.Float) 

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

    yout = scipy.zeros((len(times), len(IC)), scipy.Float)
    for ii, dyn_id in enumerate(net.dynamicVars.keys()):
        yout[:,ii] = traj.get_var_traj(dyn_id)
    pendingEvents = {}

    intervals = computeIntervals(traj)
    dv_typ_vals = [net.get_var_typical_val(dyn_id) for dyn_id
                   in net.dynamicVars.keys()]

    for start_ind, end_ind in intervals:
        start_time = times[start_ind]
        end_time = times[start_ind]
        if start_time in pendingEvents.keys():
            # We don't currently do delays in sensitivity integration, but
            #  we'll keep this here for when we do.
            IC = net.executeEvent(pendingEvents[start_time], IC, start_time)
            pendingEvents.pop(start)

        curTimes = traj.get_times()[start_ind:end_ind]

        k = min(5, end_ind - start_ind - 1)
        ys = [traj.get_var_traj(dv_id)[start_ind:end_ind] for dv_id
              in net.dynamicVars.keys()]
        tcks = [scipy.interpolate.splrep(curTimes, y, k=k, s=0) for y in ys]

        for ov_ind, opt_id in enumerate(net.optimizableVars.keys()):
            IC_this= copy.copy(IC[(ov_ind + 1) * n_dyn:
                                   (ov_ind + 2) * n_dyn])

            atol = rtol * scipy.asarray(dv_typ_vals)/\
                    net.get_var_typical_val(opt_id)

            temp = odeintr(_Ddv_and_DdvDov_dtTrunc_2,
                            IC_this, curTimes,
                            args = (net, ov_ind, tcks),
                            Dfun = Dfun_sens,
                            mxstep = 10000,
                            rtol = rtol, atol = atol)

            yout[start_ind:end_ind, (ov_ind + 1) * n_dyn:(ov_ind + 2)*n_dyn] =\
                    copy.copy(temp[0])

        IC = copy.copy(yout[end_ind - 1])

        if end_time in te:
            ind = te.index(end_time)
            event = net.events[ie[ind]]
            delay = net.fireEventAndUpdateDdv_Dov(event, yout[-1], te[ind])
            pendingEvents[te[ind] + delay] = event

    net.updateVariablesFromDynamicVars(yout[-1][:n_dyn], times[-1])

    ddv_dpTrajectory = Trajectory_mod.Trajectory(net, is_sens=True)
    ddv_dpTrajectory.appendSensFromODEINT(times, yout)

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
