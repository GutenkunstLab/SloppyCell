from SloppyCell.ReactionNetworks import *

net = Network('example')
net.add_compartment('trivial')

net.add_species('A', 'trivial', initial_conc=1)
net.add_species('B', 'trivial', initial_conc=0)
net.add_species('C', 'trivial', initial_conc=0)
net.add_species('B_scaled', 'trivial')
net.add_assignment_rule('B_scaled', 'B*1e4')

net.add_parameter('k1', 0.04)
net.add_parameter('k2', 3e7)
net.add_parameter('k3', 1e4)

net.addReaction('rxn1', stoichiometry={'A':-1,'B':+1},
                kineticLaw='k1*A')
net.addReaction('rxn2', stoichiometry={'B':-1,'C':+1},
                kineticLaw='k2*B**2')
net.addReaction('rxn3', stoichiometry={'B':-1,'A':+1},
                kineticLaw='k3*B*C')

net.set_var_optimizable('k2', False)

pred_net = net.copy()
#pred_net.set_var_ic('k2', net.get_var_ic('k2')*10)
pred_net.addReaction('rxn4', stoichiometry={'A':-1, 'C':-1, 'B':+1},
                     kineticLaw='1.0*A*C')
