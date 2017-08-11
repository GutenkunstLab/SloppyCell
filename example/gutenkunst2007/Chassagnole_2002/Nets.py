from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000051.xml', 'base')
# Remove this parameter that isn't used.
net.remove_component('time_t')
net.compile()

# We create two copies of the network, since we'll integrate it twice: once to
# focus on the long-time dynamics, and once to focus on the short time.
long = net.copy('long')
short = net.copy('short')

networks = [long, short]
int_times = [(0, 40), (0, 1)]
