import copy
import PC12Network

network = copy.deepcopy(PC12Network.baseNetwork)

network.id = 'EGFstim100'

EGFngPerml = 10002000./100.

network.setInitialVariableValue('EGF', 100*EGFngPerml)
