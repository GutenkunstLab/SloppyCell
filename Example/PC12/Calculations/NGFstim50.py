import copy
import PC12Network

network = copy.deepcopy(PC12Network.baseNetwork)

network.id = 'NGFstim50'

NGFngPerml = 228000./50.

network.setInitialVariableValue('NGF', 50*NGFngPerml)
