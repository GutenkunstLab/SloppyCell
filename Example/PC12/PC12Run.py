from SloppyCell.ReactionNetworks import *

import Nets
import Experiments as Expts
from min5_wFinal import min5_wFinal_params

m = Model([Expts.ErkMekTraverse2EGF.expt,
           Expts.Raf1LandrethEGF.expt,
           Expts.Erk1HYaoEGF.expt,
           Expts.ErkMekTraverse2NGF.expt,
           Expts.RasGreen1NGF.expt,
           Expts.Rap1YorkNGF.expt,
           Expts.Erk1HYaoNGF.expt,
           Expts.RasGreen1EGF.expt,
           Expts.Raf1BRafLandrethNGF.expt,
           Expts.ErkMekTraverseHERover.expt,
           ],
          [Nets.EGFstim100,
           Nets.NGFstim50,
           Nets.EGFstim30,
           Nets.NGFstim100,
           Nets.EGFRx50_EGFstim100,
           ])

params = KeyedList(min5_wFinal_params)
params.setOrder(m.get_params().keys())

m.compile()

c = m.cost(params)
print('Cost: {0:.5f}, should be ~39.9105'.format(c))
