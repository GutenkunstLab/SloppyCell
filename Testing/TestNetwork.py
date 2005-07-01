from SloppyCell.ReactionNetworks import *

net = Network('test')
net.add_compartment('globabl')

net.add_species('x', 'global', 1.0)
net.add_parameter('A')
net.add_rate_rule('x', '-A*x')

net.add_species('y', 'global', 2.0)
net.add_parameter('B')
net.add_rate_rule('y', '-B*y')

net.compile()
