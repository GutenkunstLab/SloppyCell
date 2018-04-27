from SloppyCell.ReactionNetworks import *

import Nets
reload(Nets)
import Experiments as Expts
import Calculations as Calcs

m = Model([Expts.Erk1HYaoEGF.expt,
           Expts.Raf1LandrethEGF.expt,
           Expts.RasGreen1EGF.expt,
           Expts.ErkMekTraverse2EGF.expt,
           Expts.Erk1HYaoNGF.expt,
           Expts.RasGreen1NGF.expt,
           Expts.ErkMekTraverse2NGF.expt,
           Expts.Rap1YorkNGF.expt,
           Expts.Raf1BRafLandrethNGF.expt,
           Expts.ErkMekTraverseHERover.expt,
           ],
          [Nets.EGFstim100,
           Nets.NGFstim50,
           Nets.EGFstim30,
           Nets.NGFstim100,
           Nets.EGFRx50_EGFstim100,
           ])

params = m.get_params().copy()
c = m.cost(params)
print('Cost before optimization: {0:.5f}, should be ~39.9105'.format(c))

perturb = 4

import numpy as np
np.random.seed(2131)
pnew = params.copy()
for ii,v in enumerate(pnew):
    pnew[ii] = v * perturb**np.random.uniform(-1,1)

print('Cost of perturbed params: {0:.5f}'.format(m.cost(pnew)))
Plotting.figure()
Plotting.plot_model_results(m, expts=['Erk1HYaoEGF', 'ErkMekTraverseHERover'])
Plotting.title('Before optimization')

# Save network with perturbed parameters to file
IO.to_SBML_file(Nets.EGFstim100, 'params_perturbed.xml')

#popt = Optimization.fmin_lm_log_params(m, pnew, disp=1, maxiter=16)

# To avoid slow optimization step during testing
popt = KeyedList([('krbEGF', 8.2921888868791371e-05), ('kruEGF', 0.027716470681392141), ('krbNGF', 1.1516917317572664e-07), ('kruNGF', 0.0012716567462045974), ('kEGF', 706.15541846106998), ('KmEGF', 6148411.044797021), ('kNGF', 439.85652370523127), ('KmNGF', 10749.326670310804), ('kdSos', 8291.9290769390936), ('KmdSos', 3169403.2607381288), ('kSos', 52.509699102554343), ('KmSos', 8031.4399723967472), ('kRasGap', 1912.9001962273949), ('KmRasGap', 820202.4948770646), ('kRasToRaf1', 0.33055186802403907), ('KmRasToRaf1', 75779.629388227127), ('kpRaf1', 389.78647896155189), ('KmpRaf1', 5792851.8278354863), ('kpBRaf', 230.74962986892515), ('KmpBRaf', 50095.148881927387), ('kdMek', 5.6528429669532558), ('KmdMek', 966755.26427205116), ('kpMekCytoplasmic', 44.730074650701155), ('KmpMekCytoplasmic', 3190639.165709862), ('kdErk', 1.9412384449826798), ('KmdErk', 604758.25742984086), ('kpP90Rsk', 0.023716435189400614), ('KmpP90Rsk', 1746126.2054943999), ('kPI3K', 29.802199423069951), ('KmPI3K', 59090.832063737049), ('kPI3KRas', 0.31598935384849713), ('KmPI3KRas', 44490.143171839023), ('kAkt', 0.020606330077855877), ('KmAkt', 142991.22649091703), ('kdRaf1ByAkt', 18.03112167994189), ('KmRaf1ByAkt', 204481.08284916033), ('kC3GNGF', 166.7296042682901), ('KmC3GNGF', 4184.2730118460477), ('kC3G', 0.75743891675182784), ('KmC3G', 6365.5804027457971), ('kRapGap', 43.733369349361638), ('KmRapGap', 874813.05813462054), ('kRap1ToBRaf', 0.59096324796478616), ('KmRap1ToBRaf', 789752.042537448), ('kdRaf1', 0.039502929899831875), ('KmdRaf1', 346.53981926951764), ('kdBRaf', 225.75098192731281), ('KmdBRaf', 5059522.129960251)])

print('Cost of optimized params: {0:.5f}'.format(m.cost(popt)))
Plotting.figure()
Plotting.plot_model_results(m, expts=['Erk1HYaoEGF', 'ErkMekTraverseHERover'])
Plotting.title('After optimization')

Plotting.show()
