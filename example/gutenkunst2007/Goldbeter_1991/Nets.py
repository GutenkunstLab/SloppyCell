from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000003.xml')
net.compile()

networks = [net]
int_times = [(0, 100)]
