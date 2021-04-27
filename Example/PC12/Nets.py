from SloppyCell.ReactionNetworks import *

# Load the network from the XML file
net = IO.from_SBML_file('BIOMD0000000033_url.xml')
net.set_var_ic('EGF', 0)
net.set_var_ic('NGF', 0)

EGFstim100 = net.copy('EGFstim100')
EGFngPerml = 100020.
EGFstim100.setInitialVariableValue('EGF', 100*EGFngPerml)

NGFstim50 = net.copy('NGFstim50')
NGFngPerml = 4560.
NGFstim50.setInitialVariableValue('NGF', 50*NGFngPerml)
