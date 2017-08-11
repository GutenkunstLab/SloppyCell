from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000015.xml')

prpp = net.copy('prpp')
prpp.set_var_ic('PRPP', net.get_var_ic('PRPP') * 10)

networks = [prpp]
int_times = [(0,70)]
