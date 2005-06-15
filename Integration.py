import copy, exceptions

import scipy

#
# Wrapper for odeint that implements support for events.
#
# func(y, t) is the rhs of the differential equations
# events is a function that returns (value, isterminal, direction) where value
#  is a list of the the condition values (event ii fires when value[ii] 
#  crosses 0). If isterminal[ii] is true the integration stops when event
#  ii fires. direction[ii] is -1, 0, or +1. If the event crossing is from + to -
#  (- to +) the event fires if direction is -1 (+1). If direction is 0 the event
#  fires whenever value[ii] crosses 0.
# eventDerivs(c, dcdt, t) is a function that returns the derivatives wrt t of
#  of the event conditions.
#
# Returns t, y, te, ye, ie. t is an array of times that includes all the 
#  requested times up to the first terminal event. y is the integration values
#  at those timepoints.
#  te, ye, ie are the times events fired at, integration values where they
#  fired, and indices of the fired events.
#
# I've tried to pass through other arguments to odeint, but I'm sure it's not
#  done well.
#
def odeintWithEvents(func, y0, t, events = None, eventDerivs = None, 
                     **kwargs):
    if kwargs.get('full_output', False):
        output, infodict = apply(scipy.integrate.odeint, (func, y0, t,), kwargs)
    else:
        output = apply(scipy.integrate.odeint, (func, y0, t,), kwargs)

    te, ye, ie = [], [], []

    # No events? Skip the rest of this rigamorole.
    if events is None or len(events(output[0], t[0])[0]) == 0:
        return t, output, te, ye, ie

    # What are the initial signs of all our events?
    value, isterminal, direction = events(output[0], t[0])
    # We round our values to the 10th decimal place in case we're restarting
    #  from when an event just fired.
    lastSigns = scipy.sign(scipy.round(value, 10))

    terminated = False
    for time, y in zip(t, output):
        value, isterminal, direction = events(y, time)
        signs = scipy.sign(scipy.round(value, 10))
        if signs != lastSigns:
            # Computes the indices of all the signs that have changed.
            changed = scipy.compress(scipy.not_equal(signs, lastSigns),
                                     range(len(signs)))
            # We don't care if a sign became non-zero, only if it changed sign
            for ii in changed:
                if lastSigns[ii] == 0:
                    changed = scipy.compress(changed != ii, changed)

            lastSigns = signs
            # What direction are our conditions changing?
            dydt = apply(func, (y, time) + kwargs.get('args', ()))
            dcdt = eventDerivs(y, dydt, time)

            #XXX: KNOWN BUG!
            #     If two terminal events fire between two time points, we'll
            #     only record the one with the lower index.
            for index in changed:
                if(direction[index] != 0
                   and direction[index] != scipy.sign(dcdt[index])):
                    continue

                # Inital conditions for our backwards integration
                y0back = scipy.concatenate((y, [0]))

                tback = scipy.array([value[index], 0])

                # This hasn't been well debugged. Basically I want to pull
                #  through any arguments that might have been passed to
                #  odeint, but the Dfun can't carry over.
                shuffledArgs = copy.copy(kwargs)
                shuffledArgs['args'] = (func, eventDerivs, index, time, 
                                       kwargs.get('args', ()))
                shuffledArgs.pop('Dfun', None)
                shuffledArgs.pop('col_deriv', None)

                # Do the backwards integration
                backoutput = apply(scipy.integrate.odeint,
                                   (funcback, y0back, tback),
                                   shuffledArgs)

                te.append(time+backoutput[-1][-1])
                ye.append(backoutput[-1][:-1])
                ie.append(index)
                if isterminal[index]:
                    terminated = True
                    break
            if terminated:
                break

    # Need to trim our outputs
    if len(te) > 0:
        t = scipy.compress(t < te[-1], t)
        output = output[:len(t)]
        t = scipy.concatenate((t, te[-1:]))
        output = scipy.concatenate((output, ye[-1:]))

    return t, output, te, ye, ie

# This is the rhs for our backwards integration.
def funcback(y, t, func, eventDerivs, index, time, args):
    dydt = apply(func, (y[:-1], time - y[-1]) + args)

    dcdt = eventDerivs(y[:-1], dydt, time - y[-1])[index]
    dydc = dydt / dcdt
    ret = scipy.concatenate((dydc, [1/dcdt]))

    return ret
