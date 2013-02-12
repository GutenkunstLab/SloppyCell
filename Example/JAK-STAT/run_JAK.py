from pylab import *                                                            # (@\label{whatever is inside these special 'at-parens comments' is visible to LaTeX}\label{code:import_start}@)
from scipy import *
from SloppyCell.ReactionNetworks import *                                      # (@\label{code:import_end}@)

net = IO.from_SBML_file('JAK-STAT_SC.xml', 'net1')                             # (@\label{code:SBML_import}@)
net.set_var_ic('v1', 'v1_0') # Won't need given initial assignments.           # (@\label{code:SBML_fix}@)

import JAK_expt                                                                # (@\label{code:import_expt}@)
m = Model([JAK_expt.expt], [net])

params = KeyedList([('r1', 0.5), ('r3', 2), ('tao', 6.0),                      # (@\label{code:params} @)
                    ('r4_0', 1.35), ('v1_0', 1.19)])

res = Residuals.PriorInLog('r3_prior', 'r3', 0, log(sqrt(1e4)))                # (@\label{code:prior_start}@)
m.AddResidual(res)
res = Residuals.PriorInLog('tao_prior', 'tao', log(4), log(sqrt(4)))
m.AddResidual(res)                                                             # (@\label{code:prior_end}@)

print 'Initial cost:', m.cost(params)                                          # (@\label{code:initial_cost}@)
params = Optimization.fmin_lm_log_params(m, params, maxiter=20, disp=False)    # (@\label{code:lm_opt}@)
print 'Optimized cost:', m.cost(params)
print 'Optimized parameters:', params

# Plot our optimal fit.
figure()                      
Plotting.plot_model_results(m)                                                 # (@\label{code:plot_results}@)

j = m.jacobian_log_params_sens(log(params))                                    # (@\label{code:j_calc}@)
jtj = dot(transpose(j), j)                                                     # (@\label{code:jtj}@)

print 'Beginning ensemble calculation.'
ens, gs, r = Ensembles.ensemble_log_params(m, asarray(params), jtj, steps=7500)# (@\label{code:ens_gen}@)
print 'Finished ensemble calculation.'

pruned_ens = asarray(ens[::25])                                                # (@\label{code:prune}@)

figure()
hist(log(pruned_ens[:,1]), normed=True)                                        # (@\label{code:hist}@)

times = linspace(0, 65, 100)
traj_set = Ensembles.ensemble_trajs(net, times, pruned_ens)                    # (@\label{code:ens_trajs}@)
lower, upper = Ensembles.traj_ensemble_quantiles(traj_set, (0.025, 0.975))

figure()
plot(times, lower.get_var_traj('frac_v3'), 'g')
plot(times, upper.get_var_traj('frac_v3'), 'g')
plot(times, lower.get_var_traj('frac_v4'), 'b')
plot(times, upper.get_var_traj('frac_v4'), 'b')

show()
