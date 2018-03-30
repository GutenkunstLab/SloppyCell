import copy
import PC12Network

network = copy.deepcopy(PC12Network.baseNetwork)

network.id = 'EGFstim30'

EGFngPerml = 3000600./30.  # Checked: same as 10002000/100 for 100 ng/ml expts

network.setInitialVariableValue('EGF', 30*EGFngPerml)
