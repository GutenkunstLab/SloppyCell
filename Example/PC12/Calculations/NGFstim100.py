import copy
import PC12Network

network = copy.deepcopy(PC12Network.baseNetwork)

network.id = 'NGFstim100'

NGFngPerml = 456000./100.  # Same as NGFngPerml = 228000./50. for 50ng/ml

network.setInitialVariableValue('NGF', 100*NGFngPerml)
