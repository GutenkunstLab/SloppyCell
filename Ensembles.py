import logging
logger = logging.getLogger('Ensembles')
import copy
import time

import scipy
import scipy.linalg
import scipy.stats

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
    f = scipy.fft.rfft(scipy.asarray(series)-scipy.mean(series), 
                       n = 2*len(series))
    # The inverse fft of |f|**2 is the autocorrelation
    ac = scipy.fft.irfft(abs(f)**2)
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
             recalc_interval=scipy.inf, recalc_func=None,
             save_hours=scipy.inf, save_to=None):
    return ensemble_log_params(m, params, hess, steps, max_run_hours,
                               temperature, step_scale, sing_val_cutoff, seeds,
                               recalc_interval, recalc_func, save_hours, 
                               save_to, log_params=False)

def ensemble_log_params(m, params, hess=None, 
                        steps=scipy.inf, max_run_hours=scipy.inf,
                        temperature=1.0, step_scale=1.0,
                        sing_val_cutoff=0, seeds=None,
                        recalc_interval=scipy.inf, recalc_func=None,
                        save_hours=scipy.inf, save_to=None, log_params=True):
    """
    Generate a Bayesian ensemble of parameter sets consistent with the data in
    the model.

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
     recalc_interval --- The sampling matrix will be recalulated each time this
                         many trial moves have been attempted.
     recalc_func --- Function used to calculate the sampling matrix. It should
                     take only a parameters argument and return the matrix.
                     If this is None, default is to use 
                     m.GetJandJtJInLogParameteters
    save_hours --- If save_to is not None, the ensemble will be saved to
                      that file every 'save_hours' hours.
    save_to --- Filename to save ensemble to.


    Outputs:
     ens, ens_costs, ratio
     ens -- List of KeyedList parameter sets in the ensemble
     ens_costs -- List of costs for each parameter set
     ratio -- Fraction of attempted moves that were accepted

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
        logger.warn('Seeding random number generator based on system time. '\
                    'Run will not be repeatable.')
    scipy.random.seed(seeds)
    if isinstance(params, KeyedList):
        param_keys = params.keys()

    curr_params = copy.deepcopy(params)
    curr_cost = m.cost(curr_params)
    ens, ens_costs = [curr_params], [curr_cost]

    # We work with arrays of params through the rest of the code
    curr_params = scipy.array(curr_params)

    if recalc_func is None:
        recalc_func = lambda p : m.GetJandJtJInLogParameters(scipy.log(p))[1]

    accepted_moves, cost_exceptions, ratio = 0, 0, scipy.nan
    start_time = last_save_time = time.time()
    samp_mat = None
    while len(ens) < steps+1:
        # Have we run too long?
        if (time.time() - start_time) >= max_run_hours*3600:
            break

        # This will always be true our first run through
        if (len(ens)%recalc_interval == 1) or (samp_mat is None):
            if (hess is None) or (len(ens) > 1):
                logger.debug('Beginning calculation of JtJ using params %s'
                             % str(ens[-1]))
                try:
                    hess = recalc_func(curr_params)
                    logger.debug('JtJ recalculated after %i steps' 
                                 % (len(ens) - 1))
                except Utility.SloppyCellException:
                    logger.warn('Calculation of new JtJ failed! ' 
                                'Continuing with previous JtJ')
            # Generate the sampling matrix used to generate candidate moves
            samp_mat, sing_vals, cutoff_sing_val = \
                    _sampling_matrix(hess, sing_val_cutoff)

        # Generate the trial move from the quadratic approximation
        deltaParams = _trial_move(samp_mat, sing_vals, cutoff_sing_val)
        # Scale the trial move by the step_scale and the temperature
        scaled_step = step_scale * scipy.sqrt(temperature) * deltaParams

        # I'm assuming log parameters here.
        if log_params:
            next_params = curr_params * scipy.exp(scaled_step)
        else:
            next_params = curr_params + scaled_step

        try:
            next_cost = m.cost(next_params)
        except Utility.SloppyCellException, X:
            logger.warn('SloppyCellException in cost evaluation at step %i, '
                        'cost set to infinity.' % len(ens))
            logger.warn('Parameters tried: %s.' % str(next_params))
            cost_exceptions += 1
            next_cost = scipy.inf

    	if _accept_move(next_cost - curr_cost, temperature):
            accepted_moves += 1.
            curr_params = next_params
            curr_cost = next_cost

        if isinstance(params, KeyedList):
            ens.append(KeyedList(zip(param_keys, curr_params)))
        else:
            ens.append(curr_params)
        ens_costs.append(curr_cost)
        ratio = accepted_moves/(len(ens) - 1)

        # Save to a file
        if save_to is not None\
           and time.time() >= last_save_time + save_hours * 3600:
            _save_ens(ens, ens_costs, ratio, save_to, cost_exceptions)
            last_save_time = time.time()

    if save_to is not None:
        _save_ens(ens, ens_costs, ratio, save_to, cost_exceptions)

    return ens, ens_costs, ratio

def _save_ens(ens, ens_costs, ratio, save_to, cost_exceptions):
    Utility.save((ens, ens_costs, ratio), save_to)
    logger.debug('Ensemble of length %i saved to %s.' % (len(ens), save_to))
    logger.debug('Acceptance ratio so far is %f.' % ratio)
    logger.debug('Cost threw an exception %i times.' % cost_exceptions)


def _accept_move(delta_cost, temperature):
    if delta_cost < 0.0:
        return True
    else:
        p = scipy.rand()
        return (p < scipy.exp(-delta_cost/temperature))

def _sampling_matrix(hessian, cutoff=0):
    ## basically need SVD of hessian - singular values and eigenvectors
    ## hessian = u * diag(singVals) * vh
    u, sing_vals, vh = scipy.linalg.svd(hessian)

    logger.debug("Hessian decomposed. Condition number is: %g."
                 % (max(sing_vals)/min(sing_vals)))

    ## scroll through the singular values and find the ones whose inverses will
    ## be huge and set them to zero also, load up the array of singular values 
    ## that we store
    ## cutoff = (1.0/_.singVals[0])*1.0e03
    ## double cutoff = _.singVals[0]*1.0e-02
    cutoff_sing_val = cutoff * max(sing_vals)

    D = 1.0/scipy.maximum(sing_vals, cutoff_sing_val)

    ## now fill in the sampling matrix ("square root" of the Hessian)
    ## note that sqrt(D[i]) is taken here whereas Kevin took sqrt(D[j])
    ## this is because vh is the transpose of his PT -JJW
    samp_mat = scipy.transpose(vh) * scipy.sqrt(D)

    return samp_mat, sing_vals, cutoff_sing_val

def _trial_move(sampling_mat, sing_vals, cutoff_sing_val):
    ## draw random numbers r ~ N(0,1)
    randVec = scipy.randn(len(sing_vals))
    # Divide the random vector by an appropriate constant such
    # that the acceptance ratio will always be about exp(-1) if
    # the harmonic approx. is good
    cutoff_vals = scipy.compress(sing_vals < cutoff_sing_val, sing_vals)
    if len(cutoff_vals):
        scale = scipy.sqrt(len(sing_vals) - len(cutoff_vals)
                           + sum(cutoff_vals)/cutoff_sing_val)
    else:
        scale = scipy.sqrt(len(sing_vals))

    randVec = randVec/scale

    ## RIGHT MULTIPLICATION
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
            traj = Dynamics.integrate(net, times, params, fill_traj=False)
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
    X = scipy.log(scipy.asarray(ens))
    X -= scipy.mean(X, 0)
    u, s, vh = scipy.linalg.svd(scipy.transpose(X))
    # This return adjust things so that can be easily compared with the JtJ
    #  eigensystem.
    return len(X)/s[::-1]**2, u[:,::-1]
