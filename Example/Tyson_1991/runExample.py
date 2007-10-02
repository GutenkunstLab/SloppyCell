import scipy
from SloppyCell.ReactionNetworks import *

# This ensures that all our plots will have smooth curves between the data 
#  points. It's a little faster to use Network.full_speed(), which only
#  calculates the points whether there is actually data.
Network.pretty_plotting()

# We access the Networks we defined in Net.py
from Nets import *

#
# First we'll reproduce figure 3
#

# Integrate all our networks.
traj1 = Dynamics.integrate(base_net, [0, 100])
traj2 = Dynamics.integrate(perturbed_net, [0,  100])
traj3 = Dynamics.integrate(growth_net, [0, 500])
traj4 = Dynamics.integrate(my_net, [0, 100])

# Figure 3 from the paper.
Plotting.figure(1, figsize=(15,4))
Plotting.subplot(1, 3, 1)
Plotting.plot_trajectory(traj1, ['M', 'YT'])
Plotting.axis([0, 100, 0, .4])
Plotting.subplot(1, 3, 2)
Plotting.plot_trajectory(traj2, ['M', 'YT'], logy=True)
Plotting.axis([0, 100, 1e-4, 1])
Plotting.subplot(2, 3, 3)
Plotting.plot_trajectory(traj3, ['k6'])
Plotting.axis([0, 500, 1, 3])
Plotting.subplot(2, 3, 6)
Plotting.plot_trajectory(traj3, ['M', 'YT'])
Plotting.axis([0, 500, 0, 0.6])

Plotting.savefig('3.png', dpi=150)

#
# Now let's fit the model to some data we've defined in Expts.py
#

from Expts import *

# We put together our Model from the experiments and networks
m = Model([expt1], [base_net, growth_net])
params = m.get_params()

print 'Initial cost:', m.cost(params)
Plotting.figure(2)
Plotting.plot_model_results(m, loc='upper right')
Plotting.title('Before fitting')

# This is a bit of python logic here...
# We loop through our parameters, and add loose priors on all their values. This
#  prevents the optimization from wandering too far.
# A prior of width log(10) corresponds to constraining the parameter to lie,
#  with 95% probility, between val/100 and val * 100.
# (These are quite tight priors. Real applications (with more constraining data)
#  will probably want looser priors.)
for id, val in params.items():
    m.AddResidual(Residuals.PriorInLog('prior_on_%s' % id, id, scipy.log(val), 
                                       scipy.log(10)))

# And now we'll optimize. Note that optimization can be tricky. One needs to
#  be in the right ball-park before the local routines can help much.
#  Once you're close, it often helps to try several of the optimizers, as
#  they have different strengths and may head different places.

# First we'll try Nelder-Mead
pmin1 = Optimization.fmin_log_params(m, params, xtol=1e-2)
# Then we run Levenburg-Marquardt
params = Optimization.leastsq_log_params(m, pmin1)
# We save our parameter values, the reload them.
Utility.save(params, 'min_params.bp')

params = Utility.load('min_params.bp')

print 'Final cost:', m.cost(params)
Plotting.figure(3)
Plotting.plot_model_results(m, loc='upper right')
Plotting.title('After fitting')

#
# Now let's run an ensemble
#

# We calculate our hessian in log parameters. 1e-4 says to vary each parameter
#  by 0.01 percent when taking the derivatives.
hess = m.hessian_log_params(params, 1e-4)

# The eigenvalues and vectors are often interesting. Let's calculate and plot
#  them.
evals, evects = Utility.eig(hess)
Plotting.figure(4)
Plotting.plot_eigvals(evals)
Plotting.figure(5)
Plotting.plot_eigvect(evects[:,0], params.keys())

# Now we'll build an ensemble of parameters.
# Make sure we run at full speed
Network.full_speed()

# Building the ensemble may generate some warnings about integration errors.
# Because our parameter space is very poorly constrained, the ensemble explores
#  regions where the equations are very difficult to integrate, so sometimes
#  the integrator won't converge. If this happens rarely, it shouldn't skew
#  the final results too much.
# In any case, the warning can be supressed with Utility.disable_warnings()
#  and re-enabled with Utility.enable_warnings().
ens, ens_costs, ratio = Ensembles.ensemble_log_params(m, params,
                                                      steps = 5000,
                                                      max_run_hours = 10/60.,
                                                      seeds=(113, 207))
Utility.save((ens, ens_costs, ratio), 'ensemble.bp')
(ens, ens_costs, ratio) = Utility.load('ensemble.bp')

# We can look at the autocorrelation function of our costs to check efficiency
Plotting.figure(7)
ac = Ensembles.autocorrelation(ens_costs)
Plotting.plot(ac)

print "Ensemble has %i members. Acceptance ratio was %f." % (len(ens), ratio)

# Let's make a prediction. How well does our data constrain the behavior of
#  growth_net?
bt, mt, st = Ensembles.net_ensemble_trajs(growth_net, 
                                          Utility.linspace(0, 500, 200),
                                          ens[::10])
Utility.save((bt, mt ,st), 'ens_trajs.bp')
best_traj, mean_traj, std_traj = Utility.load('ens_trajs.bp')
# We plot the ensemble predictions for 'CP' and 'M'.
Plotting.figure(8)
Plotting.plot_ensemble_trajs(best_traj, mean_traj, std_traj, 
                             vars=['CP', 'M'])
#
# Now for a complete data hessian
#
sens_traj = Dynamics.integrate_sensitivity(growth_net, [0, 100], fill_traj=True,
                                           params=params)
h_pd = PerfectData.hessian_log_params(sens_traj, fixed_sf=True)

# This will pop up all the plots we've generated
Plotting.show()
