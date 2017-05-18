from SloppyCell.ReactionNetworks import *

net = Network('test')
net.add_compartment('basic')

net.add_species('x', 'basic', 1.0)
net.add_parameter('A')
net.add_rate_rule('x', '-A*x')

net.add_species('y', 'basic', 2.0)
net.add_parameter('B')
net.add_rate_rule('y', '-B*y')

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

# Making a model with > 1 network so we can test parallel code
net2 = net.copy('test2')
net2.add_species('z', 'basic')
net2.add_assignment_rule('z', 'y**x + sqrt(time)')

# Make another experiment
expt2 = Experiment('expt2')
data = {'z': {0.5: (2.0, 0.2),
              2.5: (0.5, 0.05)
              }
        }
expt2.set_data({'test2': data})

m2 = Model([expt1, expt2], [net, net2])
