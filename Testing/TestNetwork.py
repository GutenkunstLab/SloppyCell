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

params = [1.0, 2.0]
    
# Make the experiment
expt1 = Experiment('expt1')
data = {'x': {1.0: (1.0, 0.1),
              2.0: (1.0, 0.1)
              },
        'y': {0.5: (2.0, 0.2),
              2.5: (0.5, 0.05)
              }
        }
expt1.set_data({'test': data})

m = Model([expt1], [net])
