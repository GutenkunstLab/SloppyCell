import scipy

from SloppyCell.ReactionNetworks import *
from SloppyCell.ReactionNetworks.RunInParallel import *

import Nets

# This Network defines our prediction
net = Nets.egf_ly
params = net.GetParameters()

N_samp = 1000        # Number of samples to use
base_Sigma = 0.5     # Sigma for well-determined parameters
missing_Sigma = 1000 # Sigma for 'missing' parameter

# Times to predict at
times = scipy.linspace(0, 45, 100)

# We build the parameter ensemble for all parameters measured
log_ens = scipy.randn(N_samp, len(params)) * scipy.log(1+0.5)/4.
log_ens += scipy.log(params)
param_ens = scipy.exp(log_ens)

# Calculate the prediction for each parameter set in the ensemble
traj_set = Ensembles.ensemble_trajs(net, times, param_ens)
# Calculate the lower and upper bounds
lower, upper = Ensembles.traj_ensemble_quantiles(traj_set, 
                                                 quantiles=(0.025, 0.975))
Utility.save((lower, upper), 'Fig.4.all.trajs.bp')
# We have to re-arrange things to plot the nice shaded regions.
xpts = scipy.concatenate((times, times[::-1]))
yy_low = lower.get_var_traj('ErkActive')/net.get_var_ic('ErkInactive')
yy_up = upper.get_var_traj('ErkActive')/net.get_var_ic('ErkInactive')
ypts_all = scipy.concatenate((yy_low, yy_up[::-1]))

# This makes the figures for all the parameters.
# To make the figure in the main text (without the collective fit bounds), use
#  for param_id in ['kpRaf1']
for param_id in params.keys():
    print param_id
    Plotting.figure(figsize=(2,1.5))
    # These parameters don't matter for this prediction
    if param_id.count('PI3K') or param_id.count('Akt') or param_id.count('C3G')\
       or param_id.count('Rap1') or param_id.count('BRaf')\
       or param_id.count('NGF'):
        continue

    log_ens = scipy.randn(N_samp, len(params)) * scipy.log(1+0.5)/4.
    log_ens += scipy.log(params)
    param_ens = scipy.exp(log_ens)

    # Generate values for the missing parameter
    missing_index = params.index_by_key(param_id)
    missing_log_ens = scipy.randn(N_samp) * scipy.log(1+missing_Sigma)/4.0
    missing_log_ens += scipy.log(params[missing_index])
    # Replace the old, tighter values that were in param_ens
    param_ens[:,missing_index] = scipy.exp(missing_log_ens)

    traj_set = Ensembles.ensemble_trajs(net, times, param_ens)
    lower, upper = Ensembles.traj_ensemble_quantiles(traj_set,
                                                 quantiles=(0.025, 0.975))

    if param_id == 'kpRaf1':
        Utility.save((lower, upper), 'Fig.4.missing.one.trajs.bp')

    yy_low = lower.get_var_traj('ErkActive')/net.get_var_ic('ErkInactive')
    yy_up = upper.get_var_traj('ErkActive')/net.get_var_ic('ErkInactive')
    ypts_miss = scipy.concatenate((yy_low, yy_up[::-1]))

    # Plot the bounds
    Plotting.fill(xpts, ypts_miss, fc='b')
    Plotting.fill(xpts, ypts_all, fc='r')

    Plotting.subplots_adjust(left=0.03, bottom=0.06, right=0.97, top=0.80)
    Plotting.axis([0, 45, 0, 1])
    Plotting.xticks([])
    Plotting.yticks([])
    Plotting.title(param_id)
    Plotting.savefig('uncert_figs/%s.eps' % param_id)
