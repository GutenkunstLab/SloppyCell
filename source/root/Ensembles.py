import logging
logger = logging.getLogger('Ensembles')
import copy
import shutil
import sys
import time

import scipy
import scipy.linalg
import scipy.stats
import scipy.fftpack

import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList
import SloppyCell.Utility as Utility

from SloppyCell import HAVE_PYPAR, my_rank, my_host, num_procs
if HAVE_PYPAR:
    import pypar

def autocorrelation(series):
    """
    Return the normalized autocorrelation of a series using the FFT.
    """
    # We need to de-mean the series. Also, we want to pad with zeros to avoid
    #  assuming our series is periodic.
    f = scipy.fftpack.rfft(scipy.asarray(series)-scipy.mean(series), 
                           n = 2*len(series))
    # The inverse fft of |f|**2 is the autocorrelation
    ac = scipy.fftpack.irfft(abs(f)**2)
    # But we padded with zeros, so it's too long
    ac = ac[:len(series)]
    # And we need to adjust to remove the effect of the zeros we added
    ac *= len(series)/(len(series) - 1.0*scipy.arange(len(series)))
    # Return the normalized ac
    return ac/ac[0]

def ensemble(m, params, hess=None, 
             steps=scipy.inf, max_run_hours=scipy.inf,
             temperature=1.0, step_scale=1.0,
             sing_val_cutoff=0, seeds=None,
             recalc_hess_alg=False, recalc_func=None,
             save_hours=scipy.inf, skip_elems=0, save_to=None,
             save_scalefactors=False):
    """
    Generate a Bayesian ensemble of parameter sets consistent with the data in
    the model. The sampling is done in terms of the bare parameters.

    Inputs:
     (All not listed here are identical to those of ensemble_log_params.)
     recalc_func --- Function used to calculate the hessian matrix. It should
                     take only a parameters argument and return the matrix.
                     If this is None, default is to use m.GetJandJtJ.
    """
    return ensemble_log_params(m, params, hess, steps, max_run_hours,
                               temperature, step_scale, sing_val_cutoff, seeds,
                               recalc_hess_alg, recalc_func, save_hours, 
                               save_to, skip_elems, log_params=False,
                               save_scalefactors=save_scalefactors)

def ensemble_log_params(m, params, hess=None, 
                        steps=scipy.inf, max_run_hours=scipy.inf,
                        temperature=1.0, step_scale=1.0,
                        sing_val_cutoff=0, seeds=None,
                        recalc_hess_alg = False, recalc_func=None,
                        save_hours=scipy.inf, save_to=None,
                        skip_elems = 0, log_params=True,
                        save_scalefactors=False):
    """
    Generate a Bayesian ensemble of parameter sets consistent with the data in
    the model. The sampling is done in terms of the logarithm of the parameters.

    Inputs:
     m -- Model to generate the ensemble for
     params -- Initial parameter KeyedList to start from 
     hess -- Hessian of the model
     steps -- Maximum number of Monte Carlo steps to attempt
     max_run_hours -- Maximum number of hours to run
     temperature -- Temperature of the ensemble
     step_scale -- Additional scale applied to each step taken. step_scale < 1
                   results in steps shorter than those dictated by the quadratic
                   approximation and may be useful if acceptance is low.
     sing_val_cutoff -- Truncate the quadratic approximation at eigenvalues
                        smaller than this fraction of the largest.
     seeds -- A tuple of two integers to seed the random number generator
     recalc_hess_alg --- If True, the Monte-Carlo is done by recalculating the
                         hessian matrix every timestep. This signficantly
                         increases the computation requirements for each step,
                         but it may be worth it if it improves convergence.
     recalc_func --- Function used to calculate the hessian matrix. It should
                     take only a log parameters argument and return the matrix.
                     If this is None, default is to use 
                     m.GetJandJtJInLogParameteters
     save_hours --- If save_to is not None, the ensemble will be saved to
                       that file every 'save_hours' hours.
     save_to --- Filename to save ensemble to.
     skip_elems --- If non-zero, skip_elems are skipped between each included 
                    step in the returned ensemble. For example, skip_elems=1
                    will return every other member. Using this option can
                    reduce memory consumption.
     save_scalefactors --- If True, scale factors will be saved during 
                          integration.

    Outputs:
     ens, ens_fes, ratio, [scale_factors]
     ens -- List of KeyedList parameter sets in the ensemble
     ens_fes -- List of free energies for each parameter set
     ratio -- Fraction of attempted moves that were accepted
     scale_factors -- List of scale factors throughout ensemble, only returned
                      if save_scalefactors is True.

    The sampling is done by Markov Chain Monte Carlo, with a Metropolis-Hasting
    update scheme. The canidate-generating density is a gaussian centered on the
    current point, with axes determined by the hessian. For a useful 
    introduction see:
     Chib and Greenberg. "Understanding the Metropolis-Hastings Algorithm" 
     _The_American_Statistician_ 49(4), 327-335
    """
    if scipy.isinf(steps) and scipy.isinf(max_run_hours):
        logger.warn('Both steps and max_run_hours are infinite! '
                    'Code will not stop by itself!')

    if seeds is None:
        seeds = int(time.time()%1 * 1e6)
        logger.debug('Seeding random number generator based on system time.')
        logger.debug('Seed used: %s' % str(seeds))
    scipy.random.seed(seeds)
    if isinstance(params, KeyedList):
        param_keys = params.keys()

    curr_params = copy.deepcopy(params)
    curr_F = m.free_energy(curr_params, temperature)
    ens, ens_Fs = [curr_params], [curr_F]
    curr_sf = m.internalVars['scaleFactors'].copy()
    ens_scale_factors = [curr_sf]

    # We work with arrays of params through the rest of the code
    curr_params = scipy.array(curr_params)

    if recalc_func is None and log_params:
        recalc_func = lambda p: m.GetJandJtJInLogParameters(scipy.log(p))[1]
    else:
        recalc_func = lambda p: m.GetJandJtJ(p)[1]

    accepted_moves, attempt_exceptions, ratio = 0, 0, scipy.nan
    start_time = last_save_time = time.time()

    # Calculate our first hessian if necessary
    if hess is None:
        hess = recalc_func(curr_params)
    # Generate the sampling matrix used to generate candidate moves
    samp_mat = _sampling_matrix(hess, sing_val_cutoff, temperature, step_scale)

    steps_attempted = 0
    while steps_attempted < steps:
        # Have we run too long?
        if (time.time() - start_time) >= max_run_hours*3600:
            break

        # Generate the trial move from the quadratic approximation
        deltaParams = _trial_move(samp_mat)
        # Scale the trial move by the step_scale and the temperature
        #scaled_step = step_scale * scipy.sqrt(temperature) * deltaParams
        scaled_step = deltaParams

        if log_params:
            next_params = curr_params * scipy.exp(scaled_step)
        else:
            next_params = curr_params + scaled_step

        try:
            next_F = m.free_energy(next_params, temperature)
        except Utility.SloppyCellException, X:
            logger.warn('SloppyCellException in free energy evaluation at step '
                        '%i, free energy set to infinity.' % len(ens))
            logger.warn('Parameters tried: %s.' % str(next_params))
            attempt_exceptions += 1
            next_F = scipy.inf
        except Utility.ConstraintViolatedException, X:
            logger.warn('ConstraintViolatedException in free energy evaluation '
                        'at step %i, free energy set to infinity.' % len(ens))
            logger.warn('Parameters tried: %s.' % str(next_params))
            attempt_exceptions += 1
            next_F = scipy.inf

        if recalc_hess_alg and not scipy.isinf(next_F):
            try:
                next_hess = recalc_func(next_params)
                next_samp_mat = _sampling_matrix(next_hess, sing_val_cutoff,
                                                 temperature, step_scale)
                accepted = _accept_move_recalc_alg(curr_F, samp_mat, 
                                                   next_F, next_samp_mat, 
                                                   deltaParams, temperature)
            except Utility.SloppyCellException, X:
                logger.warn('SloppyCellException in JtJ evaluation at step '
                            '%i, move not accepted.' % len(ens))
                logger.warn('Parameters tried: %s.' % str(next_params))
                attempt_exceptions += 1
                next_F = scipy.inf
                accepted = False
        else:
            accepted = _accept_move(next_F - curr_F, temperature)

        steps_attempted += 1
        if accepted:
            accepted_moves += 1.
            curr_params = next_params
            curr_sf = m.internalVars['scaleFactors'].copy()
            curr_F = next_F
            if recalc_hess_alg:
                hess = next_hess
                samp_mat = next_samp_mat

        if steps_attempted % (skip_elems + 1) == 0:
            ens_Fs.append(curr_F)
            if save_scalefactors:
                ens_scale_factors.append(curr_sf)
            if isinstance(params, KeyedList):
                ens.append(KeyedList(zip(param_keys, curr_params)))
            else:
                ens.append(curr_params)
        ratio = accepted_moves/steps_attempted

        # Save to a file
        if save_to is not None\
           and time.time() >= last_save_time + save_hours * 3600:
            _save_ens(ens, ens_Fs, ratio, save_to, attempt_exceptions,
                      steps_attempted, ens_scale_factors, 
                      save_sf=save_scalefactors)
            last_save_time = time.time()

    if save_to is not None:
        _save_ens(ens, ens_Fs, ratio, save_to, attempt_exceptions, 
                  steps_attempted, ens_scale_factors,
                  save_sf=save_scalefactors)

    if save_scalefactors:
        return ens, ens_Fs, ratio, ens_scale_factors
    else:
        return ens, ens_Fs, ratio

def _save_ens(ens, ens_Fs, ratio, save_to, attempt_exceptions, steps_attempted,
              ens_scale_factors, save_sf=False):
    temp_name = save_to + '_temporary_SloppyCellFile'
    if not save_sf:
        Utility.save((ens, ens_Fs, ratio), temp_name)
    else:
        Utility.save((ens, ens_Fs, ratio, ens_scale_factors), temp_name)
    shutil.move(temp_name, save_to)
    logger.debug('Ensemble of length %i saved to %s.' % (len(ens), save_to))
    logger.debug('Acceptance ratio so far is %f.' % ratio)
    logger.debug('Attempted moves threw an exception %i times out of %i ' 
                 'attempts.' % (attempt_exceptions, steps_attempted))

def _accept_move(delta_F, temperature):
    """
    Basic Metropolis accept/reject step.
    """
    p = scipy.rand()
    return (p < scipy.exp(-delta_F/temperature))

def _accept_move_recalc_alg(curr_F, curr_samp_mat, next_F, next_samp_mat, 
                            step, T):
    """
    Accept/reject when each the sampling matrix is recalculated each step.
    """
    pi_x = scipy.exp(-curr_F/T)
    # This is the current location's covariance sampling matrix
    sigma_curr = scipy.dot(curr_samp_mat, scipy.transpose(curr_samp_mat))
    sigma_curr_inv = scipy.linalg.inv(sigma_curr)
    # This is the transition probability from the current point to the next.
    q_x_to_y = scipy.exp(-_quadratic_cost(step, sigma_curr_inv))\
            / scipy.sqrt(scipy.linalg.det(sigma_curr))

    pi_y = scipy.exp(-next_F/T)
    sigma_next = scipy.dot(next_samp_mat, scipy.transpose(next_samp_mat))
    sigma_next_inv = scipy.linalg.inv(sigma_next)
    q_y_to_x = scipy.exp(-_quadratic_cost(-step, sigma_next_inv))\
            / scipy.sqrt(scipy.linalg.det(sigma_next))

    p = scipy.rand()
    accepted = (pi_y*q_y_to_x)/(pi_x*q_x_to_y)
    did_accepted = p<accepted

    return p < accepted


def _sampling_matrix(hessian, cutoff=0, temperature=1, step_scale=1):
    # basically need SVD of hessian - singular values and eigenvectors
    # hessian = u * diag(singVals) * vh
    u, sing_vals, vh = scipy.linalg.svd(0.5 * hessian)

    # scroll through the singular values and find the ones whose inverses will
    # be huge and set them to zero also, load up the array of singular values 
    # that we store
    # cutoff = (1.0/_.singVals[0])*1.0e03
    # double cutoff = _.singVals[0]*1.0e-02
    cutoff_sing_val = cutoff * max(sing_vals)

    D = 1.0/scipy.maximum(sing_vals, cutoff_sing_val)

    ## now fill in the sampling matrix ("square root" of the Hessian)
    ## note that sqrt(D[i]) is taken here whereas Kevin took sqrt(D[j])
    ## this is because vh is the transpose of his PT -JJW
    samp_mat = scipy.transpose(vh) * scipy.sqrt(D)

    # Divide the sampling matrix by an additional factor such
    # that the expected quadratic increase in cost will be about 1.
    cutoff_vals = scipy.compress(sing_vals < cutoff_sing_val, sing_vals)
    if len(cutoff_vals):
        scale = scipy.sqrt(len(sing_vals) - len(cutoff_vals)
                           + sum(cutoff_vals)/cutoff_sing_val)
    else:
        scale = scipy.sqrt(len(sing_vals))

    samp_mat /= scale
    samp_mat *= step_scale
    samp_mat *= scipy.sqrt(temperature)

    return samp_mat

def _trial_move(sampling_mat):
    randVec = scipy.randn(len(sampling_mat))
    trialMove = scipy.dot(sampling_mat, randVec)

    return trialMove

def _quadratic_cost(trialMove, hessian):
    """
    The cost from the quadratic approximation of a trialMove, given the hessian.

    (Note: the hessian here is assumed to be the second derivative matrix of the
     cost, without an additional factor of 1/2.)
    """
    quadratic = 0.5*scipy.dot(scipy.transpose(trialMove), 
                              scipy.dot(hessian, trialMove))
    return quadratic

def traj_ensemble_stats(traj_set):
    """
    Return the mean and standard deviation trajectory objects for the given 
    input list of trajectories. (All must be evaluated at the same points.)
    """
    all_values = [traj.values for traj in traj_set]

    mean_values = scipy.mean(all_values, 0)
    std_values = scipy.std(all_values, 0)

    mean_traj = copy.deepcopy(traj_set[0])
    mean_traj.values = mean_values

    std_traj = copy.deepcopy(traj_set[0])
    std_traj.values = std_values

    return mean_traj, std_traj

def few_ensemble_trajs(net, times, elements):
    import SloppyCell.ReactionNetworks.Dynamics as Dynamics
    traj_set = []
    for params in elements:
        try:
            traj = Dynamics.integrate(net, times, params=params, 
                                      fill_traj=False)
            if not scipy.any(scipy.isnan(traj.values)):
                traj_set.append(traj)
        except Utility.SloppyCellException:
            logger.warn('Exception in network integration on node %i.'
                        % my_rank)

    return traj_set

def ensemble_trajs(net, times, ensemble):
    """
    Return a list of trajectories evaluated at times for all parameter sets
    in ensemble.
    """
    traj_set = []
    elems_assigned = [ensemble[node::num_procs] for node in range(num_procs)]
    for worker in range(1, num_procs):
        command = 'Ensembles.few_ensemble_trajs(net, times, elements)'
        args = {'net': net, 'times': times, 'elements': elems_assigned[worker]}
        pypar.send((command, args), worker)

    traj_set = few_ensemble_trajs(net, times, elems_assigned[0])

    for worker in range(1, num_procs):
        traj_set.extend(pypar.receive(worker))

    return traj_set

def net_ensemble_trajs(net, times, ensemble):
    traj_set = ensemble_trajs(net, times, ensemble)
    best_traj = traj_set[0]
    mean_traj, std_traj = traj_ensemble_stats(traj_set)
    return best_traj, mean_traj, std_traj

def traj_ensemble_quantiles(traj_set, quantiles=(0.025, 0.5, 0.975)):
    """
    Return a list of trajectories, each one corresponding the a given passed-in
    quantile.
    """
    all_values = scipy.array([traj.values for traj in traj_set])
    sorted_values = scipy.sort(all_values, 0)
                   
    q_trajs = []
    for q in quantiles:
        # Calculate the index corresponding to this quantile. The q is because
        #  Python arrays are 0 indexed
        index = q * (len(sorted_values) - 1)
        below = int(scipy.floor(index))
        above = int(scipy.ceil(index))
        if above == below:
            q_values = sorted_values[below]
        else:
            # Linearly interpolate...
            q_below = (1.0*below)/(len(sorted_values)-1)
            q_above = (1.0*above)/(len(sorted_values)-1)
            q_values = sorted_values[below] + (q - q_below)*(sorted_values[above] - sorted_values[below])/(q_above - q_below)
        q_traj = copy.deepcopy(traj_set[0])
        q_traj.values = q_values
        q_trajs.append(q_traj)

    return q_trajs

def PCA_eig_log_params(ens):
    """
    Return the Principle Component Analysis eigenvalues and eigenvectors (in 
     log parameters) of an ensemble. (This function takes the logs for you.)
    """
    return PCA_eig(scipy.log(scipy.asarray(ens)))

def PCA_eig(ens):
    """
    Return the Principle Component Analysis eigenvalues and eigenvectors 
     of an ensemble.
    """
    X = scipy.asarray(ens)
    mean_vals = scipy.mean(X, 0)
    X -= mean_vals
    u, s, vh = scipy.linalg.svd(scipy.transpose(X))
    # This return adjusts things so that can be easily compared with the JtJ
    #  eigensystem.
    X += mean_vals
    return len(X)/s[::-1]**2, u[:,::-1]
