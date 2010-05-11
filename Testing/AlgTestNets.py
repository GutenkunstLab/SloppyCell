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
# The reaction 'out' has the stoichiometryMath expression 1/1
# for X1 to test the stoichiometry math compatibility of SloppyCell
algebraic_net.addReaction('out', kineticLaw = 'k2 * S2',
                          stoichiometry = {'X1': '1/1', 'T': -1,'S2':0.0})

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

# algebraic_net_multi is a system that has 3 algebraic rules
# the dynamics should be the same as in algebraic_net, but
# having extra algebraic rules helps test whether systems with many algebraic
# equations will work correctly.

algebraic_net_multi = algebraic_net.copy('alg_net_multi')
algebraic_net_multi.add_parameter('S3', initial_value = 0.2, is_constant = False)
algebraic_net_multi.add_parameter('S4', initial_value = 0.2, is_constant = False)
algebraic_net_multi.add_algebraic_rule('S3-k1')
algebraic_net_multi.add_algebraic_rule('S4-k2')

# algebraic_net_andor_events has trigger and event assignments that make use
# of and_ and or_funcs with multiple arguments
# also uses a variable set to 'time' to make sure that SBML time output is OK.

algebraic_net_andor_events = algebraic_net.copy('algebraic_net_andor_events')
algebraic_net_andor_events.add_parameter('timevar', is_constant=False)
eas = KeyedList()
piecewise_X0 = 'piecewise(0, or_func(gt(X0, 100), gt(X1, 100), lt(T, 0)), 1)'
piecewise_X1 = 'piecewise(1, or_func(gt(X0, 100), gt(X1, 100), lt(T, 0)), 0)'
eas.set('X0', piecewise_X0)
eas.set('X1', piecewise_X1)
eas.set('timevar', 'time')
logical_trigger = 'and_func(gt(1, 0), lt(1, 0), eq(10, 11))'
algebraic_net_andor_events.addEvent(id='logical_event', trigger=logical_trigger,
                                    eventAssignments=eas)

# algebraic_net_under has an extra unconstrained variable that appears
# in an algebraic rule, so SloppyCell should complain that the system
# is underconstrained.

algebraic_net_under = algebraic_net.copy('alg_net_underdetermined')
algebraic_net_under.addParameter('P', value=5.0, isConstant=False)
del algebraic_net_under.algebraicRules[0]
algebraic_net_under.add_algebraic_rule('S2 + S1 - T*(P/P)')

# algebraic_net_constraints implements a test of the Constraint
# objects.  Also add some events to make sure constraints don't
# interfere with events.

algebraic_net_constraints = algebraic_net.copy('alg_net_constraints')

eventAssignments=KeyedList()
eventAssignments.set('X0',0.1)
algebraic_net_constraints.addEvent('event1', 'and(gt(time,18),gt(X1,0.4))',eventAssignments)
algebraic_net_constraints.addEvent('event2', 'and(gt(time,10),gt(X1,0.4))',eventAssignments)
algebraic_net_constraints.addConstraint('X1toobig',trigger='leq(X1,0.5)',
                                        message='X1 is big!')

# algebraic_net_manual is the same as algebraic net, except we make the cross
# referencing manual before creating the network. _makeCrossReferences(...) must
# be called before the network is integrated.

algebraic_net_manual = Network('algebraic_net_manual')
algebraic_net_manual._manualCrossReferences(flag=True)

algebraic_net_manual.add_compartment('cell')
algebraic_net_manual.add_species('X0', 'cell')
algebraic_net_manual.add_species('X1', 'cell')
algebraic_net_manual.add_species('T', 'cell')
algebraic_net_manual.add_species('S1', 'cell')
algebraic_net_manual.add_species('S2', 'cell')

algebraic_net_manual.add_parameter('Keq')
algebraic_net_manual.add_parameter('k1')
algebraic_net_manual.add_parameter('k2')

algebraic_net_manual.add_algebraic_rule('S2 + S1 - T')
algebraic_net_manual.add_assignment_rule('S2', 'Keq * S1')

algebraic_net_manual.addReaction('in', kineticLaw = 'k1 * X0',
                          stoichiometry = {'X0': -1, 'T': 1})
# The reaction 'out' has the stoichiometryMath expression 1/1
# for X1 to test the stoichiometry math compatibility of SloppyCell
algebraic_net_manual.addReaction('out', kineticLaw = 'k2 * S2',
                          stoichiometry = {'X1': '1/1', 'T': -1,'S2':0.0})

algebraic_net_manual.set_var_ics({'X0': 1.0,
                           'X1': 0.0,
                           'T': 0.0,
                           'S1': 0.0,
                           'S2': 0.0,
                           'Keq': 2.5,
                           'k1': 0.1,
                           'k2': 0.15})

