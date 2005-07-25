import copy, sets, sys, types
import scipy
import SloppyCell.Integration as Integration

# need these just for the situation that an initial condition for
# the sensitivity equations depends on some messy function of parameters
# Maybe the initial condition should be evaluated in Network and passed in?
log, log10 = scipy.log, scipy.log10
exp = scipy.exp
cos, sin, tan = scipy.cos, scipy.sin, scipy.tan
acos, asin, atan = scipy.arccos, scipy.arcsin, scipy.arctan
cosh, sinh, tanh = scipy.cosh, scipy.sinh, scipy.tanh
arccosh, arcsinh, arctanh = scipy.arccosh, scipy.arcsinh, scipy.arctanh
exponentiale, pi = scipy.e, scipy.pi

def _getIntervals(net, timepoints):
    cT = sets.Set([timepoints[0], timepoints[-1]] +
                  [e.triggeringTime for e in net.timeTriggeredEvents])
    cT = scipy.sort(list(cT))
    cT = scipy.compress(cT >= timepoints[0] and cT <= timepoints[-1], cT)
    return zip(cT[:-1], cT[1:])

def Integrate(net, timepoints, rtol = None):
    timepoints = scipy.array(timepoints)

    # We want to stop the integrator whenever we have a timeTriggered event.
    intervals = _getIntervals(net, timepoints)

    output = scipy.zeros((0, len(net.dynamicVars)), scipy.Float)
    outTimes, outTE, outYE, outIE = [], [], [], []
    pendingEvents = {}

    IC = net.getDynamicVarValues()
    func = net.get_ddv_dt
    Dfun = net.get_d2dv_ddvdt

    atol = None
    if scipy.isscalar(rtol):
        atol = [rtol*var.typicalValue for var in net.dynamicVars.values()]
    elif rtol is not None:
        rtol = None
        print >> sys.stderr, "Non-scalar rtol passed into Integrate(). Haven't decided how to handle this yet, so we're defaulting to rtol = None, atol = None (odeint defaults)."

    if getattr(net, 'integrateWithLogs', False):
        def func(logDV, time):
            dv = scipy.exp(logDV)
            ddv_dt = net.get_ddv_dt(dv, time)
            return ddv_dt / dv
        Dfun = None
        atol = 1e-12
        rtol = 1e-6
        
    for (start, intervalEnd) in intervals:
        # Check whether any timeTriggered events need to fire
        for id, event in net.timeTriggeredEvents.items():
            if event.triggeringTime == start:
                outIE.append(net.events.indexByKey(id))
                outYE.append(IC)
                outTE.append(start)
                delay = net.fireEvent(event, IC, start)
                pendingEvents[round(start + delay, 8)] = event

        while start < intervalEnd:
            # Executing any pending events
            if round(start, 8) in pendingEvents.keys():
                IC = net.executeEvent(pendingEvents[round(start, 8)], IC, start)
                pendingEvents.pop(round(start, 8))

            # What's the farthest our integration can possibly run?
            end = min(pendingEvents.keys() + [intervalEnd])
            curPoints = scipy.compress(scipy.logical_and(timepoints > start,
                                                         timepoints < end),
                                       timepoints)
            curPoints = scipy.concatenate(([start], curPoints, [end]))


            if getattr(net, 'integrateWithLogs', False):
                IC = scipy.log(IC)
            # Do the integration, getting initial conditions from the net
            t, odeint_array, te, ye, ie = Integration.\
                    odeintWithEvents(func, copy.copy(IC), curPoints, 
                                     events = net.get_eventValues,
                                     eventDerivs = net.get_eventDerivs,
                                     Dfun = Dfun, col_deriv = 1,
                                     mxstep = 10000, rtol = rtol, atol = atol)

            if getattr(net, 'integrateWithLogs', False):
                output = scipy.concatenate((output, scipy.exp(odeint_array)))
            else:
                output = scipy.concatenate((output, odeint_array))

            outTimes.extend(t)

            start, IC = outTimes[-1], copy.copy(output[-1])

            outTE.extend(te)
            outYE.extend(ye)
            outIE.extend(ie)

            # If a complex event fired, add it to the pending list
            if len(te) > 0:
                event = net.complexEvents[ie[-1]]
                delay = net.fireEvent(event, odeint_array[-1], te[-1])
                pendingEvents[round(te[-1] + delay, 8)] = event

    # Update the network with the final variable values
    net.updateVariablesFromDynamicVars(output[-1], t)
    return outTimes, output, outTE, outYE, outIE

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


def Integrate_Ddv_Dov(net, timepoints, rtol = None):
    #XXX: Need to handle case where initial values are functions of parameters
    timepoints = scipy.array(timepoints)

    # We want to stop the integrator whenever we have a timeTriggered event.
    intervals = _getIntervals(net, timepoints)

    nDyn, nOpt = len(net.dynamicVars), len(net.optimizableVars)

    atolDv = None
    if scipy.isscalar(rtol):
        atolDv = scipy.array([rtol*var.typicalValue for
                              var in net.dynamicVars.values()])
    elif rtol is not None:
        rtol = None
        print >> sys.stderr, "Non-scalar rtol passed into Integrate_Ddv_Dov(). Haven't decided how to handle this yet, so we're defaulting to rtol = None, atol = None (odeint defaults)."

    output = scipy.zeros((0, nDyn * (nOpt + 1)), scipy.Float)
    pendingEvents = {}

    IC = scipy.zeros(nDyn * (nOpt + 1), scipy.Float) # IC is zero unless IC is a function of parameters
    dvInd = 0
    for id, var in net.dynamicVars.items():
	if type(var.initialValue) == types.StringType:
		rule = net.substituteFunctionDefinitions(var.initialValue)
		for ovInd,ovName in enumerate(net.optimizableVars.keys()) :
			DwrtOV = net.takeDerivative(rule, ovName)
			if DwrtOV != '0' :
				DwrtOV = net.substituteFunctionDefinitions(DwrtOV)
				DwrtOV = net.substituteVariableNames(DwrtOV)
				DwrtOV = DwrtOV.replace('self','net')
				IC[nDyn*(ovInd+1) + dvInd] = eval(DwrtOV)
	dvInd = dvInd + 1

    IC[:nDyn] = net.getDynamicVarValues()
    outTimes, outTE, outYE, outIE = [], [], [], []

    for (start, intervalEnd) in intervals:
        # Check whether any timeTriggered events need to fire
        for id, event in net.timeTriggeredEvents.items():
            if event.triggeringTime == start:
                outIE.append(net.events.indexByKey(id))
                outYE.append(IC)
                outTE.append(start)
                delay = net.fireEventAndUpdateDdv_Dov(event, IC, start)
                pendingEvents[round(start + delay, 8)] = event

        while start < intervalEnd:
            # Executing any pending events
            if round(start, 8) in pendingEvents.keys():
                IC = net.executeEventAndUpdateDdv_Dov(pendingEvents
                                                      [round(start, 8)], 
                                                      IC, start)
                pendingEvents.pop(round(start, 8))

            # What's the farthest our integration can possibly run?
            end = min(pendingEvents.keys() + [intervalEnd])
            curPoints = scipy.compress(scipy.logical_and(timepoints > start,
                                                         timepoints < end),
                                       timepoints)
            curPoints = scipy.concatenate(([start], curPoints, [end]))

            # We'll do an initial integration run without the sensitivites
            #  to find out how far we need to integrate to get to a terminating
            #  event. This reduces the amount of backtracking we need to do
            #  when doing the much slower sensitivity integrations.
            t, odeint_array, te, ye, ie = Integration.\
                    odeintWithEvents(net.get_ddv_dt, copy.copy(IC[:nDyn]),
                                     curPoints, 
                                     events = net.get_eventValues,
                                     eventDerivs = net.get_eventDerivs,
                                     Dfun=net.get_d2dv_ddvdt, col_deriv=1,
                                     mxstep=10000, rtol = rtol, atol = atolDv)

            newSpace = scipy.zeros((len(t), nDyn * (nOpt + 1)), scipy.Float)
            output = scipy.concatenate((output, newSpace)) 
            output[-len(t):, :nDyn] = copy.copy(odeint_array[:, :nDyn])
            outTimes.extend(t)

            for ovIndex in range(nOpt):
                ICTrunc = scipy.zeros(2 * nDyn, scipy.Float)
                ICTrunc[:nDyn] = copy.copy(IC[:nDyn])
                ICTrunc[nDyn:] = copy.copy(IC[(ovIndex + 1) * nDyn:
                                              (ovIndex + 2) * nDyn])

                # Absolute tolerances for sensitivity variables
                if atolDv is not None:
                    atolSens = atolDv/net.optimizableVars[ovIndex].typicalValue
		    atolForThis = scipy.concatenate((atolDv, atolSens))
                else:
                    atolForThis = None

                odeint_array = scipy.integrate.odeint(Ddv_and_DdvDov_dtTrunc, 
                                                      ICTrunc, t, 
                                                      args = (net, ovIndex),
                                                      mxstep = 10000, 
                                                      rtol = rtol, 
                                                      atol = atolForThis)

                output[-len(t):, (ovIndex + 1) * nDyn:(ovIndex + 2)*nDyn] =\
                        copy.copy(odeint_array[:, nDyn:])

            start, IC = outTimes[-1], output[-1]

            outTE.extend(te)
            outYE.extend(ye)
            outIE.extend(ie)

            if len(te) > 0:
                event = net.complexEvents[ie[-1]]
                delay = net.fireEventAndUpdateDdv_Dov(event, output[-1], te[-1])
                pendingEvents[round(te[-1] + delay, 8)] = event

    net.updateVariablesFromDynamicVars(output[-1][:nDyn], t)
    return outTimes, output, outTE, outYE, outIE
