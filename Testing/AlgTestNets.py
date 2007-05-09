"""
Implementation of networks in algebraicRules-assignment_in_algebraic.xml and
algebraicRules-fastReactionExample-l2.xml from the SBML semantic test suite.
"""
from SloppyCell.ReactionNetworks import *

algebraic_net = Network('algebraic_net')

algebraic_net.add_compartment('cell')
algebraic_net.add_species('X0', 'cell')
algebraic_net.add_species('X1', 'cell')
algebraic_net.add_species('T', 'cell')
algebraic_net.add_species('S1', 'cell')
algebraic_net.add_species('S2', 'cell')

algebraic_net.add_parameter('Keq')
algebraic_net.add_parameter('k1')
algebraic_net.add_parameter('k2')

algebraic_net.add_algebraic_rule('S2 + S1 - T')
algebraic_net.add_assignment_rule('S2', 'Keq * S1')

algebraic_net.addReaction('in', kineticLaw = 'k1 * X0',
                          stoichiometry = {'X0': -1, 'T': 1})
algebraic_net.addReaction('out', kineticLaw = 'k2 * S2',
                          stoichiometry = {'X1': 1, 'T': -1})

algebraic_net.set_var_ics({'X0': 1.0,
                           'X1': 0.0,
                           'T': 0.0,
                           'S1': 0.0,
                           'S2': 0.0,
                           'Keq': 2.5,
                           'k1': 0.1,
                           'k2': 0.15})

algebraic_net_assignment = algebraic_net.copy('alg_net_assigned')
algebraic_net_assignment.add_parameter('S_sum')
algebraic_net_assignment.add_parameter('Alg_Rule_RHS')
algebraic_net_assignment.add_parameter('p_1', is_constant=False)
algebraic_net_assignment.add_parameter('p_2', is_constant=False)
algebraic_net_assignment.add_parameter('p_3', is_constant=False)

algebraic_net_assignment.add_assignment_rule('S_sum', 'S1 + S2')
algebraic_net_assignment.add_assignment_rule('Alg_Rule_RHS', 'S_sum - T')
del algebraic_net_assignment.algebraicRules[0]
algebraic_net_assignment.add_algebraic_rule('Alg_Rule_RHS')

algebraic_net_assignment.add_event('event0', 'lt(X0, 0.5)', {'p_2': 1})
algebraic_net_assignment.set_var_ics({'p_1': 0,
                                      'p_2': 0,
                                      'p_3': 0})
