from SloppyCell.ReactionNetworks import *

# Here we have an experiment where data was taken on M in two conditions, 
#  corresponding to the networks 'base', and 'perturbed'
expt1 = Experiment('expt1')
expt1.SetData({'base':{
                       'M':{10.0: (0.2, 0.1),
                            20.0: (0.05, 0.03),
                            45.0: (0.2, 0.1)},
                       'CP':{45.0: (1.2, 0.2),
                             5.0: (0.9, 0.1),
                             35.0: (0.85, 0.05),
                             }
                       },
               'growth':{
                         'M':{20.0: (0.02, 0.03),
                              87.0: (0.2, 0.1),
                              150.0: (0.03, 0.02),
                              200.0: (0.2, 0.1)}
                         }
               })
