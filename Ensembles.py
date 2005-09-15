import copy
import scipy
import time

import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

eig = scipy.linalg.eig

set_seeds = scipy.stats.seed
set_seeds(72529486,916423761)

def autocorrelation(series):
    """
    Return the autocorrelation of a series
    """
    series = scipy.array(series)
    slen = len(series)
    smean = scipy.mean(series)
    acorr = scipy.zeros(len(series), scipy.Float)
    c0 = scipy.sum((series - smean)**2)/slen

    acorr[0] = 1.0
    for ii in range(1, len(acorr)):
        acorr[ii] = scipy.sum((series[:-ii]-smean) * (series[ii:]-smean))\
                /(slen*c0)

    return acorr

def ensemble_log_params(m, params, hess, 
                        steps=scipy.inf, max_run_hours=scipy.inf,
                        temperature=1.0, step_scale=1.0,
                        sing_val_cutoff=1e-4, seeds=None):
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
        raise ValueError, 'Both steps and max_run_hours cannot be infinity, or the code will never stop!'

    if seeds is not None:
        set_seeds(seeds[0], seeds[1])

    curr_params = params.copy()
    param_keys = curr_params.keys()
    curr_cost = m.cost(curr_params)
    ens, ens_costs = [curr_params], [curr_cost]

    curr_params = scipy.array(curr_params)

    # Generate the sampling matrix used to generate candidate moves
    samp_mat, sing_vals, cutoff_sing_val = _sampling_matrix(hess, 
                                                            sing_val_cutoff)

    accepted_moves = 0
    start_time = time.time()

    while len(ens) < steps+1:
        # Have we run too long?
        if (time.time() - start_time)/3600 > max_run_hours:
            break

        # Generate the trial move from the quadratic approximation
        deltaParams = _trial_move(samp_mat, sing_vals, cutoff_sing_val)
        # Scale the trial move by the step_scale and the temperature
        scaled_step = step_scale * scipy.sqrt(temperature) * deltaParams

        # I'm assuming log parameters here.
        next_params = curr_params * scipy.exp(scaled_step)
        # If we had non-log parameters, it would be:
        # next_params = curr_params + scaled_step

        next_cost = m.cost(next_params)
    	if _accept_move(next_cost - curr_cost, temperature):
            accepted_moves += 1.
            curr_params = next_params
            curr_cost = next_cost

        ens.append(KeyedList(zip(param_keys, curr_params)))
        ens_costs.append(curr_cost)

    ratio = accepted_moves/len(ens)
    return ens, ens_costs, ratio

def _accept_move(delta_cost, temperature):
    if delta_cost < 0.0:
        return True
    else:
        p = scipy.rand()
        return (p < scipy.exp(-delta_cost/temperature))

def _sampling_matrix(hessian, cutoff=1e-4):
    ## basically need SVD of hessian - singular values and eigenvectors
    ## hessian = u * diag(singVals) * vh
    u, sing_vals, vh = scipy.linalg.svd(hessian)

    print "Hessian decomposed. Condition number is: %g "  % \
            (max(sing_vals)/min(sing_vals))

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
    quadratic = scipy.dot(scipy.transpose(trialMove), 
                          scipy.dot(hessian, trialMove))
    return quadratic

def net_ensemble_trajs(net, times, ensemble):
    best_traj = net.integrate(times, ensemble[0], addTimes=False)

    all_trajs = [best_traj]
    for params in ensemble:
        all_trajs.append(net.integrate(times, params, addTimes=False))

    all_values = [traj.values for traj in all_trajs]

    mean_values = scipy.mean(all_values, 0)
    std_values = scipy.std(all_values, 0)

    mean_traj = copy.deepcopy(best_traj)
    mean_traj.values = mean_values

    std_traj = copy.deepcopy(best_traj)
    std_traj.values = std_values

    return best_traj, mean_traj, std_traj
