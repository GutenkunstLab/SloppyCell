import copy
import sets

import scipy

from SloppyCell.lsodar import odeintr
from SloppyCell.ReactionNetworks import KeyedList
import Trajectory

def integrate(net, times, params=None, rtol=1e-6):
    if params is not None:
        net.update_optimizable_vars(params)

    if getattr(net, 'times_to_add', None):
        times = sets.Set(times)
        times.union_update(sets.Set(net.times_to_add))
        times = list(times)
        times.sort()

    times = scipy.array(times)
    if net.add_tail_times:
        times = scipy.concatenate((times, [1.05*times[-1]]))

    # If you ask for time = 0, we'll assume you want dynamic variable values
    #  reset.
    if times[0] == 0:
        net.resetDynamicVariables()

    # te, ye, and ie are the times, dynamic variable values, and indices
    #  of fired events
    IC = net.getDynamicVarValues()
    start = times[0]
    func = net.get_ddv_dt

    atol = None
    if scipy.isscalar(rtol):
        atol = [rtol*var.typicalValue for var in net.dynamicVars.values()]
    elif rtol is not None:
        if len(rtol) == len(IC):
            atol = [rtol[ii] *var.typicalValue for ii, var 
                    in enumerate(net.dynamicVars.values())]
        else:
            sys.stderr.write("Non-scalar rtol of improper length passed into integrate!")

    if getattr(net, 'integrateWithLogs', False):
        def func(logDV, time):
            dv = scipy.exp(logDV)
            ddv_dt = net.get_ddv_dt(dv, time)
            return ddv_dt / dv
        atol = 1e-12
        rtol = 1e-6

    yout, tout = scipy.zeros((0, len(IC)), scipy.Float), []
    #yout, tout = scipy.array([IC]), [start]
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

        nRoots = len(net.events)
        if getattr(net, 'integrateWithLogs', False):
            Dfun_temp = None
            IC = scipy.log(IC)
        else:
            Dfun_temp = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))
        temp = odeintr(func, copy.copy(IC), curTimes, 
                       root_func = net.root_func,
                       root_term = [True]*nRoots,
                       Dfun = Dfun_temp, 
                       mxstep = 10000, rtol = rtol, atol = atol,
                       int_pts = net.add_int_times,
                       insert_events = True,
                       full_output = True)

        if getattr(net, 'integrateWithLogs', False):
            yout = scipy.concatenate((yout, scipy.exp(temp[0])))
        else:
            yout = scipy.concatenate((yout, temp[0]))
        tout.extend(temp[1])
        start, IC = tout[-1], copy.copy(yout[-1])
        te.extend(temp[2])
        ye.extend(temp[3])
        ie.extend(temp[4])

        # If an event fired
        if temp[4]:
            event = net.events[ie[-1]]
            delay = net.fireEvent(event, yout[-1], te[-1])
            pendingEvents[te[-1] + delay] = event

    net.updateVariablesFromDynamicVars(yout[-1], tout[-1])

    if not net.add_int_times:
        yout = _reduce_times(yout, tout, times)
        tout = times

    indexKeys = net.dynamicVars.keys() + net.assignedVars.keys()
    keyToColumn = KeyedList(zip(indexKeys, range(len(indexKeys))))

    trajectory = Trajectory.Trajectory(net, keyToColumn)
    trajectory.appendFromODEINT(tout, yout)

    return trajectory

def integrate_sensitivity(net, times, params=None, rtol=1e-6):
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
        atolDv = [rtol*var.typicalValue for var in net.dynamicVars.values()]
    elif rtol is not None:
        if len(rtol) == len(IC):
            atolDv = [rtol[ii] *var.typicalValue for ii, var 
                    in enumerate(net.dynamicVars.values())]
        else:
            sys.stderr.write("Non-scalar rtol of improper length passed into integrate_sensitivity!")

    start = times[0]
    # IC is zero unless IC is a function of parameters
    IC = scipy.zeros(nDyn * (nOpt + 1), scipy.Float) 

    # Handle sensitivities of initial conditions
    for dvInd, (id, var) in enumerate(net.dynamicVars.items()):
	if isinstance(var.initialValue, str):
            rule = net.substituteFunctionDefinitions(var.initialValue)
            for ovInd, ovName in enumerate(net.optimizableVars.keys()) :
                DwrtOV = net.takeDerivative(rule, ovName)
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

        for ovIndex in range(nOpt):
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
                atolSens = scipy.asarray(atolDv)\
                        /net.optimizableVars[ovIndex].typicalValue
	        atolForThis = scipy.concatenate((atolDv, atolSens))
            else:
                atolForThis = None

            temp2 = odeintr(Ddv_and_DdvDov_dtTrunc,
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

    alldvs = net.dynamicVars.keys() + net.assignedVars.keys()
    keyToColumnSensNames = alldvs + [(cname, pname) for
                                     pname in net.optimizableVars.keys()
                                     for cname in alldvs]
    keyToColumnSens = [(name, index) for index, name in
                       enumerate(keyToColumnSensNames)]
    keyToColumnSens = KeyedList(keyToColumnSens)
    
    ddv_dpTrajectory = Trajectory.Trajectory(net, keyToColumnSens)
    ddv_dpTrajectory.appendSensFromODEINT(tout, yout)

    return ddv_dpTrajectory

def D2dv_Dov_dtTrunc(y, time, net, ovIndex):
    nDV = len(net.dynamicVars)

    # The dynamic variables
    dv = y[:nDV]
    # The current total derivates of the dynamic vars wrt the optimizable var
    #  whose index is ovIndex
    Ddv_Dov = y[nDV:]

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

def Ddv_and_DdvDov_dtTrunc(y, time, net, ovIndex):
    nDV = len(net.dynamicVars)

    # XXX: Might want to fill in a global array here. Perhaps a small speed
    #      boost
    Dc_Dt = net.get_ddv_dt(y[:nDV], time)
    D2dv_Dov_dt = D2dv_Dov_dtTrunc(y, time, net, ovIndex)

    return scipy.concatenate((Dc_Dt, D2dv_Dov_dt))

def _reduce_times(yout, tout, times):
    jj = 0
    for ii, twanted in enumerate(times):
        while tout[jj] != twanted:
            jj += 1
        yout[ii] = yout[jj]

    yout = yout[:len(times)]

    return yout
