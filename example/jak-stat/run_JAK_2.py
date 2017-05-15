from pylab import *                                                            # (@\label{whatever is inside these special 'at-parens comments' is visible to LaTeX}\label{code:import_start}@)
from scipy import *
from SloppyCell.ReactionNetworks import *                                      # (@\label{code:import_end}@)

net = IO.from_SBML_file('JAK-STAT_SC.xml', 'net1')                             # (@\label{code:SBML_import}@)
net.set_var_ic('v1', 'v1_0') # Won't need given initial assignments.           # (@\label{code:SBML_fix}@)

import JAK_expt2                                                                # (@\label{code:import_expt}@)
m = Model([JAK_expt2.expt], [net])

params = KeyedList([('r1', 0.1), ('r3', 3), ('tao', 6.0),                      # (@\label{code:params} @)
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
show()
