import copy
import PC12Network

network = copy.deepcopy(PC12Network.baseNetwork)

network.id = 'EGFRx50_EGFstim100'

EGFngPerml = 10002000./100.

network.setInitialVariableValue('EGF', 100*EGFngPerml)
network.setInitialVariableValue('freeEGFReceptor', 
                                50*PC12Network.baseNetwork.getInitialVariableValue('freeEGFReceptor'))
