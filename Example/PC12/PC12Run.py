from SloppyCell.ReactionNetworks import *

import Nets
import Experiments as Expts

m = Model([
           Expts.ErkMekTraverse2EGF.expt,
           Expts.ErkMekTraverse2NGF.expt,
           Expts.Raf1LandrethEGF.expt,
           Expts.Rap1YorkNGF.expt,
           Expts.RasGreen1NGF.expt,
           ],
          [Nets.EGFstim100,
           Nets.NGFstim50,
           ])

params = m.get_params().copy()
c = m.cost(params)
print('Cost before optimization: {0:.5f}'.format(c))

perturb = 4

import numpy as np
np.random.seed(2131)
pnew = params.copy()
for ii,v in enumerate(pnew):
    pnew[ii] = v * perturb**np.random.uniform(-1,1)

print('Cost of perturbed params: {0:.5f}'.format(m.cost(pnew)))
Plotting.figure()
Plotting.plot_model_results(m)
Plotting.title('Before optimization')

## Save network with perturbed parameters to file
#IO.to_SBML_file(Nets.EGFstim100, 'params_perturbed.xml')

popt = Optimization.fmin_lm_log_params(m, pnew, disp=1, maxiter=8)

print('Cost of optimized params: {0:.5f}'.format(m.cost(popt)))
Plotting.figure()
Plotting.plot_model_results(m)
Plotting.title('After optimization')

Plotting.show()
