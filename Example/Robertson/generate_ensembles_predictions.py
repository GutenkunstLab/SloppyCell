from numpy import *
import scipy.stats
from SloppyCell.ReactionNetworks import *
import example_net, example_model

#
# Calculate an ensemble
#
# Step scale = 3 gives faster convergence for this problem.
print("Generating ensemble")
ens_data, costs, ratio = Ensembles.ensemble_log_params(example_model.m, 
                                                       example_model.popt, 
                                                       hess=example_model.jtj, 
                                                       steps=10000, 
                                                       step_scale=3)
ens_data = array(ens_data)
Utility.save(ens_data, 'example.ens_data.bpkl')

#
# Calculate a cost surface
# Note that we're defining our cost here to *not* include terms from priors,
# so we create a copy of our model without the priors and use that.
#
print("Generating cost surface")
Npt = 41
land_xx = logspace(-4, 2, Npt)
land_yy = logspace(1.5, 7.5, Npt)
Z = zeros((len(land_yy), len(land_xx)))
for ii,y in enumerate(land_yy):
    for jj,x in enumerate(land_xx):
        p = [x,y]
        Z[ii,jj] = example_model.noprior_m.cost(p)
Utility.save((land_xx, land_yy, Z), 'example.model_surface.bpkl')

ens_data = Utility.load('example.ens_data.bpkl')

#
# Sample from "all measured" ensemble
#
measure_uncert = 0.25
Nsamp = len(ens_data)
ens_all_measured = exp(log(example_model.popt) + measure_uncert*random.randn(Nsamp,2))

#
# Sample from "one unmeasured" ensemble
#
ens_one_measured = zeros((Nsamp, 2))
ens_one_measured[:,0] = exp(log(example_model.popt[0]) + measure_uncert*random.randn(Nsamp))
ens_one_measured[:,1] = exp(log(example_model.popt[1]) + example_model.jtj_uncerts[1]*random.randn(Nsamp))

#
# Sample from jtj hessian
#
samp_mat = Ensembles._sampling_matrix(example_model.jtj)
ens_jtj = array([Ensembles._trial_move(samp_mat) for ii in range(Nsamp)])
ens_jtj = exp(log(example_model.popt) + ens_jtj)

Utility.save((ens_data, ens_all_measured, ens_one_measured, ens_jtj), 
             'example.ensembles.bpkl')

tt = linspace(0,50,100)
def calc_trajectories(ens, tt):
    yys = []
    for p in ens:
        traj = Dynamics.integrate(example_net.pred_net, tt, params=p, 
                                  fill_traj=False)
        yy = traj.get_var_traj('C')
        yys.append(yy)
    return yys

def calc_trajectory_bounds(yys, lw_percent=2.5, up_percent=97.5):
    up, lw = [], []
    for yt in array(yys).transpose():
        up.append(scipy.stats.scoreatpercentile(yt, 97.5))
        lw.append(scipy.stats.scoreatpercentile(yt, 2.5))
    return array(lw), array(up)

#
# Calculate predictions
#
pred_data = calc_trajectories(ens_data[::10], tt)
pred_all_measured = calc_trajectories(ens_all_measured[::10], tt)
pred_one_measured = calc_trajectories(ens_one_measured[::10], tt)
pred_jtj = calc_trajectories(ens_jtj[::10], tt)

#
# Calculate 95% bounds on predictions
#
pred_data_lw, pred_data_up = calc_trajectory_bounds(pred_data)
pred_all_measured_lw, pred_all_measured_up = calc_trajectory_bounds(pred_all_measured)
pred_one_measured_lw, pred_one_measured_up = calc_trajectory_bounds(pred_one_measured)
pred_jtj_lw, pred_jtj_up = calc_trajectory_bounds(pred_jtj)

Utility.save((pred_data_lw, pred_data_up, 
              pred_all_measured_lw, pred_all_measured_up,
              pred_one_measured_lw, pred_one_measured_up, 
              pred_jtj_lw, pred_jtj_up), 'example.pred_bounds.bpkl')
