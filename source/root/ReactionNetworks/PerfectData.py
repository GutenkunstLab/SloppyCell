"""
Methods for generating "perfect data" hessians.
"""
__docformat__ = "restructuredtext en"

import copy
import logging
logger = logging.getLogger('RxnNets.PerfectData')

import sets

import scipy
import scipy.integrate

import Dynamics

from SloppyCell import HAVE_PYPAR, my_rank, my_host, num_procs
if HAVE_PYPAR:
    import pypar

def apply_func_to_traj(traj, func, only_nonderivs=False):
    """
    Return a trajectory with func applied to each variable stored in the
    trajectory
    """
    if only_nonderivs:
        keys = [key for key in self.key_column.keys()
                if not isinstance(key, tuple)]
    else:
        keys = None

    ret_traj = traj.copy_subset(keys)
    for var, col in traj.key_column.items():
        vals = func(traj, var)
        ret_traj.values[:,col] = vals

    return ret_traj 

def update_typical_vals(networks, int_times, rtol = 1e-9, fraction=1.0,
                        cutoff=1e-14): 
    """
    Update the typical var values for a group of networks.

    Find the maximum of each variable over the integrations. In each network
    the typical value is set to fraction of that maximum. If that maximum value
    is less than cutoff, the typical value is set to 1.

    networks    List of networks to work with
    int_times   List of corresponding integration endpoints 
                  (ie. [(0, 100), (0, 50)])
    fraction    Relative size of typical value, compared to max value over
                the integrations.
    rtol        Relative tolerance for the integrations.
    cutoff      Values below this are considered to be zero
    """
    max_vals = {}
    for net, times in zip(networks, int_times):
        traj = Dynamics.integrate(net, times, rtol=rtol, fill_traj=True)
        for var_id in net.variables.keys():
            curr_max = max(traj.get_var_traj(var_id))
            max_vals[var_id] = max(curr_max, max_vals.get(var_id, 0))

    for var_id, val in max_vals.items():
        for net in networks:
            if net.variables.has_key(var_id):
                if val > cutoff:
                    net.set_var_typical_val(var_id, val*fraction)
                else:
                    net.set_var_typical_val(var_id, 1.0)

    return max_vals

def typ_val_uncert(fraction = 0.1, cutoff=1e-14):
    """
    This is an uncertainty that is fraction of the variable's typical value.
    """
    def sigmaFunc(traj, data_id):
        sigma = traj.get_var_typical_val(data_id) * fraction
        if sigma < cutoff:
            logger.warn('sigma < cutoff value (%g) for variable %s! '
                        'Taking sigma = 1.' % (cutoff, data_id))
            sigma = 1
        return sigma
    return sigmaFunc

def discrete_data(net, params, pts, interval, vars=None, random=False,
                  uncert_func=typ_val_uncert(0.1, 1e-14)):
    """
    Return a set of data points for the given network generated at the given
    parameters.

    net         Network to generate data for
    params      Parameters for this evaluation of the network
    pts         Number of data points to output
    interval    Integration interval
    vars        Variables to output data for, defaults to all species in net
    random      If False data points are distributed evenly over interval
                If True they are spread randomly and uniformly over each
                variable
    uncert_func Function that takes in a trajectory and a variable id and
                returns what uncertainty should be assumed for that variable,
                either as a scalar or a list the same length as the trajectory.
    """
    # Default for vars
    if vars is None:
        vars = net.species.keys()

    # Assign observed times to each variable
    var_times = {}
    for var in vars:
        if random:
            var_times[var] = scipy.rand(pts) * (interval[1]-interval[0]) + interval[0]
        else:
            var_times[var] = scipy.linspace(interval[0], interval[1], pts)

    # Create a sorted list of the unique times in the var_times dict
    int_times = sets.Set(scipy.ravel(var_times.values()))
    int_times.add(0)
    int_times = list(int_times)
    int_times.sort()

    # Get the trajectory
    traj = Dynamics.integrate(net, int_times, params=params, fill_traj=False)

    # Build up our data dictionary
    data = {}
    for var, times in var_times.items():
        var_data = {}
        data[var] = var_data

        # Calculate our uncertainties
        var_uncerts = uncert_func(traj, var)
        for time in times:
            val = traj.get_var_val(var, time)
            if scipy.isscalar(var_uncerts):
                uncert = var_uncerts
            else:
                index = traj._get_time_index(time)
                uncert = var_uncerts[index]
            var_data[time] = (val, uncert)
    return data

def hessian_log_params(sens_traj, data_ids=None, opt_ids=None, 
                       fixed_sf=False, return_dict=False, 
                       uncert_func=typ_val_uncert(1.0, 1e-14)):
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
    uncert_func Function that takes in a trajectory and a variable id and
                returns what uncertainty should be assumed for that variable,
                either as a scalar or a list the same length as the trajectory.
    """
    if data_ids is None:
        data_ids = sens_traj.dynamicVarKeys + sens_traj.assignedVarKeys
    if opt_ids is None:
        opt_ids = sens_traj.optimizableVarKeys

    data_sigmas = {}
    for data_id in data_ids:
        ds = uncert_func(sens_traj, data_id)
        if scipy.isscalar(ds):
            ds = scipy.zeros(len(sens_traj), scipy.float_) + ds
        data_sigmas[data_id] = ds

    vars_assigned = [data_ids[node::num_procs] for node in range(num_procs)]
    for worker in range(1, num_procs):
        logger.debug('Sending to worker %i.' % worker)
        # reduce the amount we have to pickle
        # The only things the worker needs in the sens_traj are those that
        #  refer to data_ids it has to deal with.
        vars_needed = sets.Set(sens_traj.optimizableVarKeys)
        vars_needed.union_update(vars_assigned[worker])
        for var in vars_assigned[worker]:
            vars_needed.union_update([(var, ov) for ov in opt_ids])
        worker_traj = sens_traj.copy_subset(vars_needed)
        # And the only uncertainties it needs have to do with those data_ids
        worker_ds = dict([(var, data_sigmas[var])
                          for var in vars_assigned[worker]])
        command = 'PerfectData.compute_sf_LMHessian_conts(sens_traj, data_ids,'\
                'data_sigmas, opt_ids, fixed_sf)'
        args = {'sens_traj': worker_traj, 'data_ids': vars_assigned[worker], 
                'data_sigmas': worker_ds, 'opt_ids': opt_ids,
                'fixed_sf': fixed_sf}
        pypar.send((command, args), worker)

    hess_dict = compute_sf_LMHessian_conts(sens_traj, vars_assigned[0],
                                           data_sigmas, opt_ids, fixed_sf)

    for worker in range(1, num_procs):
        logger.debug('Receiving from worker %i.' % worker)
        hess_dict.update(pypar.receive(worker))

    hess = scipy.sum(hess_dict.values(), axis=0)
    if return_dict:
        return hess, hess_dict
    else:
        return hess

def compute_sf_LMHessian_conts(sens_traj, data_ids, data_sigmas, opt_ids,
                               fixed_sf):
    hess_dict = {}
    for data_id in data_ids:
        if fixed_sf:
            sf_deriv = dict([(id, 0) for id in opt_ids])
        else:
            sf_deriv = get_sf_derivs(sens_traj, data_id, data_sigmas[data_id],
                                     opt_ids)
        hess_dict[data_id] = computeLMHessianContribution(sens_traj, data_id, 
                                                          data_sigmas[data_id], 
                                                          opt_ids, sf_deriv)

    return hess_dict

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
        if intTheorySq != 0:
            scaleFactorDerivs[optId] /= intTheorySq
        else:
            scaleFactorDerivs[optId] = 0

    return scaleFactorDerivs

def computeLMHessianContribution(traj, dataId, data_sigma, optIds, 
                                 scaleFactorDerivs):
    LMHessian = scipy.zeros((len(optIds), len(optIds)), scipy.float_)

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
