import copy
import sys

import scipy

import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

def fmin_powell_log_params(m, params, *args, **kwargs):
    func = m.cost_log_params

    pmin = scipy.optimize.fmin_powell(func, scipy.log(params), 
                                      *args, **kwargs)
    if isinstance(params, KeyedList):
        pout = params.copy()
        pout.update(scipy.exp(pmin))
        return pout
    else:
        return scipy.exp(pmin)

def fmin_log_params(m, params, *args, **kwargs):
    func = m.cost_log_params

    pmin = scipy.optimize.fmin(func, scipy.log(params),
                               *args, **kwargs)

    if isinstance(params, KeyedList):
        pout = params.copy()
        pout.update(scipy.exp(pmin))
        return pout
    else:
        return scipy.exp(pmin)

def leastsq_log_params(m, params, *args, **kwargs):
    func = m.res_log_params

    pmin, msg = scipy.optimize.leastsq(func, scipy.log(params), *args, **kwargs)

    if isinstance(params, KeyedList):
        pout = params.copy()
        pout.update(scipy.exp(pmin))
        return pout
    else:
        return scipy.exp(pmin)
