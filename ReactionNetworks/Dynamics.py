import copy

import scipy

import SloppyCell.lsodar as lsodar
from SloppyCell.ReactionNetworks import KeyedList, Trajectory

def integrate(net, times, params = None):
    if params is not None:
        net.update_optimizable_vars(params)

    times = scipy.array(times)
    if net.natural_times:
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

    #yout = scipy.zeros((0, len(net.dynamicVars)), scipy.Float)
    yout, tout = scipy.array([IC]), [start]
    te, ye, ie = [], [], []
    pendingEvents = {}

    rtol, atol = None, None

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
        Dfun_temp = lambda y, t: scipy.transpose(net.get_d2dv_ddvdt(y, t))
        temp = lsodar.odeintr(func, copy.copy(IC), curTimes, 
                              root_func = net.root_func,
                              root_term = [True]*nRoots,
                              Dfun = Dfun_temp, 
                              mxstep = 10000, rtol = rtol, atol = atol,
                              int_pts = net.natural_times,
                              insert_events = net.natural_times,
                              full_output = True)

        IC = copy.copy(temp[0][-1])
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

    if not net.natural_times:
        jj = 0
        for ii, twanted in enumerate(times):
            while tout[jj] != twanted:
                jj += 1
            yout[ii] = yout[jj]

        yout = yout[:len(times)]
        tout = times

    indexKeys = net.dynamicVars.keys() + net.assignedVars.keys()
    keyToColumn = KeyedList(zip(indexKeys, range(len(indexKeys))))

    trajectory = Trajectory.Trajectory(net, keyToColumn)
    trajectory.appendFromODEINT(tout, yout)

    return trajectory
