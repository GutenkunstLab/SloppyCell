"""
Implementation of network in stoichiometry-mathML-l2.xml from
the SBML semantic test suite.
"""
from SloppyCell.ReactionNetworks import *

stoichMath_net = Network('stoichMath_net')

stoichMath_net.add_compartment('c1')
stoichMath_net.add_species('A', 'c1')
stoichMath_net.add_species('B', 'c1')

stoichMath_net.add_parameter('k')

stoichMath_net.addReaction('r1', kineticLaw = 'A * k',
                          stoichiometry = {'A': '-B*10', 'B': 1})

stoichMath_net.set_var_ics({'A': 1.0,
                           'B': 0.0,
                           'k': 1.0})

