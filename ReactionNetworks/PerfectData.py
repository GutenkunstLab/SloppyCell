"""
Methods for generating "perfect data" hessians.
"""
__docformat__ = "restructuredtext en"

import logging
logger = logging.getLogger('RxnNets.PerfectData')

import scipy

_SIGMA_CUTOFF = 1e-14

# This specifies the uncertainties assumed.
def sigmaFunc(traj, data_id):
    """
    This takes the uncertainty in a variable to be equal to its typical value 
    divided by 10.
    """
    sigma = traj.get_var_typical_val(data_id)/10.0
    if sigma < _SIGMA_CUTOFF:
        logger.warn('sigma < cutoff value (%g) for variable %s! '
                    'Taking sigma = 1.' % (_SIGMA_CUTOFF, data_id))
        sigma = 1

    return sigma

def hessian_log_params(sens_traj, data_ids=None, opt_ids=None, fixed_sf=False,
                       return_dict=False):
    """
    Calculate the "perfect data" hessian in log parameters given a sensitivity
    trajectory.

    sens_traj   Sensitivity trajectory of Network being considered.
    data_ids    A sequence of variable id's to assume we have data for. If 
                data_ids is None, all dynamic and assigned variables will be 
                used.
    opt_ids     A sequence of parameter id's to calculate derivatives with 
                respect to. The hessian is (len(opt_ids) x len(opt_ids)).
                If opt_ids is None, all optimizable variables are considered.
    fixed_sf    If True, calculate the hessian assuming fixed scale factors.
    return_dict If True, returned values are (hess, hess_dict). hess_dict is a
                dictionary keyed on the elements of data_ids; each corresponding
                value is the hessian assuming data only on a single variable.
                hess is the sum of all these hessians
    """
    if data_ids is None:
        data_ids = sens_traj.dynamicVarKeys + sens_traj.assignedVarKeys
    if opt_ids is None:
        opt_ids = sens_traj.optimizableVarKeys

    sf_derivs, hess_dict = {}, {}
    for data_id in data_ids:
        data_sigma = sigmaFunc(sens_traj, data_id)
        if scipy.isscalar(data_sigma):
            data_sigma = scipy.zeros(len(sens_traj), scipy.Float) + data_sigma

        if fixed_sf:
            sf_derivs[data_id] = dict([(id, 0) for id in opt_ids])
        else:
            sf_derivs[data_id] = get_sf_derivs(sens_traj, data_id, data_sigma, 
                                               opt_ids)

        hess_dict[data_id] = \
                computeLMHessianContribution(sens_traj, data_id, data_sigma, 
                                             opt_ids, sf_derivs[data_id])

    hess = scipy.sum(hess_dict.values())
    if return_dict:
        return hess, hess_dict
    else:
        return hess

def get_intervals(traj):
    # We want to break up our integrals when events fire, so first we figure out
    #  when they fired by looking for duplicated times in the trajectory
    times = traj.get_times()
    eventIndices = scipy.compress(scipy.diff(times) == 0, 
                                  scipy.arange(len(times)))
    intervals = zip([0] + list(eventIndices + 1), 
                    list(eventIndices + 1) + [len(times)])

    return intervals

def get_sf_derivs(traj, dataId, data_sigma, optIds):
    scaleFactorDerivs = {}
    intTheorySq = 0
    for start, end in get_intervals(traj):
        y = traj.get_var_traj(dataId)[start:end]
        times = traj.get_times()[start:end]
        sigma = data_sigma[start:end]

        value = scipy.integrate.simps((y/sigma)**2, times, 
                                      even='last')
        intTheorySq += value
    
        for optId in optIds:
            optValue = abs(traj.get_var_traj(optId)[start:end])
            sens = traj.get_var_traj((dataId, optId))[start:end]*optValue

            numerator = scipy.integrate.simps(sens*y/sigma**2, times, 
                                              even='last')

            scaleFactorDerivs.setdefault(optId, 0)
            scaleFactorDerivs[optId] += -numerator

    for optId in optIds:
        scaleFactorDerivs[optId] /= intTheorySq

    return scaleFactorDerivs

def computeLMHessianContribution(traj, dataId, data_sigma, optIds, 
                                 scaleFactorDerivs):
    LMHessian = scipy.zeros((len(optIds), len(optIds)), scipy.Float)

    # We break up our integral at event firings.
    for start, end in get_intervals(traj):
        times = traj.timepoints[start:end]
        y = traj.getVariableTrajectory(dataId)[start:end]
        sigma = data_sigma[start:end]

        for optIndex1, optId1 in enumerate(optIds):
            # We convert our sensitivity trajectory to a sensitivity wrt the 
            #  log(abs()) by multiplying by the parmeter value.
            optValue = abs(traj.get_var_traj(optId1)[start:end])
            sens1 = traj.get_var_traj((dataId, optId1))[start:end]*optValue
            dB1 = scaleFactorDerivs[optId1]
            for jj, optId2 in enumerate(optIds[optIndex1:]):
                optIndex2 = jj + optIndex1

                optValue = abs(traj.get_var_traj(optId2)[start:end])
                sens2 = traj.get_var_traj((dataId, optId2))[start:end]*optValue
                dB2 = scaleFactorDerivs[optId2]

                integrand = (sens1 + dB1 * y) * (sens2 + dB2 * y)/ sigma**2
                # We do the even='last' for speed. Otherwise, for an even number
                #  of points, simps does twice as much work as for an odd
                #  number.
                # In tests it really doesn't make much difference in accurary.
                value = scipy.integrate.simps(integrand, times, 
                                              even='last')

                LMHessian[optIndex1][optIndex2] += value
                if optIndex1 != optIndex2:
                    LMHessian[optIndex2][optIndex1] += value
    
    LMHessian /= (traj.timepoints[-1] - traj.timepoints[0])

    return LMHessian
