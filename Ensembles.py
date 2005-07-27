import copy
import scipy

from SloppyCell import KeyedList

eig = scipy.linalg.eig

set_seeds = scipy.stats.seed
set_seeds(72529486,916423761)

def build_ensemble_log_params(m, params, hess, ens_size, 
                              step_scale = 1.0, temperature = 1.0):
    """
    Build an ensemble of parameter sets for the given model.
    """
    samp_mat, sing_vals, cutoff_sing_val = _sampling_matrix(hess)

    param_keys = params.keys()

    curr_params = scipy.array(params)
    orig_cost = m.cost(curr_params)
    ens_params, ens_costs = [copy.deepcopy(params)], [orig_cost]
    ens_count, trial_moves = 1, 0
    while ens_count < ens_size:
        deltaParams = _trial_move(samp_mat, sing_vals, cutoff_sing_val)
        quadCost = _quadratic_cost(deltaParams, hess)

        # I'm assuming log parameters here.
        next_params = curr_params * scipy.exp(step_scale * deltaParams)
        # If we had non-log parameters, it would be:
        # next_params = curr_params + step_scal*deltaParams

        next_cost = m.cost(next_params)
    	trial_moves += 1
    	if _accept_move(quadCost, next_cost - orig_cost, temperature):
            ens_count += 1
            ens_params.append(KeyedList(zip(param_keys, next_params)))
            ens_costs.append(next_cost)
            curr_params = next_params

    print 'Acceptance ratio: %f' % ((ens_count - 1.0)/trial_moves)
    return ens_params, ens_costs

def _accept_move(quadCost, deltaCost, temp):
    if deltaCost < 0.0:
        return True
    else:
        p = scipy.rand()
        return (p < scipy.exp(-deltaCost/temp))

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
    scale = scipy.sqrt(len(sing_vals) - len(cutoff_vals)
                       + sum(cutoff_vals)/cutoff_sing_val)
    randVec = randVec/scale
    ## now rescale by the sampling matrix, including the appropriate
    ## reweighting by the control parameter.

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

    lower_traj, upper_traj = copy.deepcopy(best_traj), copy.deepcopy(best_traj)
    lower_traj.values = mean_values - std_values
    upper_traj.values = mean_values + std_values

    return best_traj, mean_traj, lower_traj, upper_traj


def GetEnsembleTrajs(calc, times):
    # pass in networkcalculation (an instance of the network)
    # , and a vector of times. 
    # Get out a dictionary indexed by chemical
    # which has the envelope of all the trajectories of all
    # the chemicals over the ensemble. Also passes out the
    # variance at each time point so can compare with the
    # linearized prediction for the variance
    
    bestfittraj = {}
    meanchem = {}
    variancechems = {}
    maxtrajchems = {}
    mintrajchems = {}
    varstocalc = {}
    alldvs = calc.dynamicVars.keys() + calc.assignedVars.keys()  # all chemicals
    for var in alldvs :
    	maxtrajchems[var] = scipy.zeros((len(times),),scipy.Float)
    	mintrajchems[var] = scipy.ones((len(times),),scipy.Float)*1.0e80
    	meanchem[var] = scipy.zeros((len(times),),scipy.Float)
    	variancechems[var] = scipy.zeros((len(times),),scipy.Float)
    	bestfittraj[var] = []
    	varstocalc[var] = times
    
    
    for var in alldvs :
    	# store best fit separately
    	for params in _.ensembleParams[0:1] :
    		calc.Calculate(varstocalc,params)
    		tmptraj = calc.trajectory.getVariableTrajectory(var)
    		bestfittraj[var] = scipy.array(tmptraj)
    
    	# this computes trajectories for all parameter sets in the ensemble
    	# May want to step through this
    	for params in _.ensembleParams :
    		calc.Calculate(varstocalc,params)
    		tmptraj = calc.trajectory.getVariableTrajectory(var)
    		# alltrajs[var].append(tmptraj)
    		meanchem[var] = meanchem[var] + scipy.array(tmptraj)
    		variancechems[var] = variancechems[var] + scipy.array(tmptraj)**2 # an elementwise square
    		for i in range(0,len(times)) :
    			maxtrajchems[var][i] = max(maxtrajchems[var][i],tmptraj[i])
    			mintrajchems[var][i] = min(mintrajchems[var][i],tmptraj[i])
    
    for kys in meanchem.keys() :
    	meanchem[kys] = meanchem[kys]/len(_.ensembleParams)
    	variancechems[kys] = variancechems[kys]/len(_.ensembleParams) - meanchem[kys]**2
    
    return bestfittraj,meanchem,variancechems,maxtrajchems,mintrajchems
