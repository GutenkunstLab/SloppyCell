import copy
import sys

import scipy

def fmin_powell(m, params, print_costs=False, *args, **kwargs):
    if print_costs:
        def func(log_params):
            c = m.cost_log_params(log_params)
            print 'Cost: %g' % c
            sys.stdout.flush()
            return c
    else:
        func = m.cost_log_params

    pmin = scipy.optimize.fmin_powell(func, scipy.log(params), *args, **kwargs)
    pout = copy.deepcopy(params)
    pout.update(scipy.exp(pmin))
    return pout

def fmin(m, params, print_costs=False, *args, **kwargs):
    if print_costs:
        def func(log_params):
            c = m.cost_log_params(log_params)
            print 'Cost: %g' % c
            sys.stdout.flush()
            return c
    else:
        func = m.cost_log_params

    pmin = scipy.optimize.fmin(func, scipy.log(params), *args, **kwargs)
    pout = copy.deepcopy(params)
    pout.update(scipy.exp(pmin))
    return pout

def leastsq(m, params, *args, **kwargs):
    func = m.res_log_params

    pmin, msg = scipy.optimize.leastsq(func, scipy.log(params), *args, **kwargs)
    pout = copy.deepcopy(params)
    pout.update(scipy.exp(pmin))
    return pout
