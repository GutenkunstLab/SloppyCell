from SloppyCell.ReactionNetworks import *

# Load the network from the XML file
net = IO.from_SBML_file('BIOMD0000000033.xml')
net.set_var_ic('EGF', 0)
net.set_var_ic('NGF', 0)
net.compile()

EGFstim100 = net.copy('EGFstim100')
EGFngPerml = 100020.
EGFstim100.setInitialVariableValue('EGF', 100*EGFngPerml)

EGFstim30 = net.copy('EGFstim30')
EGFstim30.setInitialVariableValue('EGF', 30*EGFngPerml)

NGFstim50 = net.copy('NGFstim50')
NGFngPerml = 4560.
NGFstim50.setInitialVariableValue('NGF', 50*NGFngPerml)

NGFstim100 = net.copy('NGFstim100')
NGFstim100.setInitialVariableValue('NGF', 100*NGFngPerml)

EGFRx50_EGFstim100 = net.copy('EGFRx50_EGFstim100')
EGFRx50_EGFstim100.setInitialVariableValue('EGF', 100*EGFngPerml)
EGFRx50_EGFstim100.setInitialVariableValue('freeEGFReceptor',
                                           50*net.getInitialVariableValue('freeEGFReceptor'))
