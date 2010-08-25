import unittest
import os
import copy

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *
# Check whether we actually have the SBML methods.
_HAVE_SBML = (hasattr(IO, 'to_SBML_file') and hasattr(IO, 'from_SBML_file'))

# Load the fast reaction example from the SBML semantic test suite.
# To avoid extra dependencies on libsbml, we use verions built by SloppyCell.
from AlgTestNets import algebraic_net, algebraic_net_assignment, algebraic_net_multi,\
     algebraic_net_andor_events
tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 51)])
sbml_file = os.path.join('SBML_files','algebraicRules-fastReactionExample-l2.xml')

class test_SBMLInterface(unittest.TestCase):
    def test_tofile(self):
        """ Basic test of SBML output """

        # First remove the output file in case it was produced in a
        # previous test
        if os.path.exists('algebraic_net.xml') == True:
            os.remove('algebraic_net.xml')

        # Now output the SBML and make sure it got created        
        SBMLInterface.toSBMLFile(algebraic_net, 'algebraic_net.xml')

        self.assertEqual(os.path.exists('algebraic_net.xml'), True)

    def test_fromfile(self):
        """ Basic test of SBML input """

        fastReaction_net = SBMLInterface.fromSBMLFile(sbml_file)
        
        self.assertEqual(fastReaction_net.id, 'algebraicRules_fastReactionExample')

    def test_sbml_with_time(self):
        """ Test of SBML output/input using time in an event assignment """

        # First remove the output file in case it was produced in a
        # previous test
        if os.path.exists('algebraic_net_time.xml') == True:
            os.remove('algebraic_net_time.xml')

        # Now output the SBML and make sure it got created        
        SBMLInterface.toSBMLFile(algebraic_net_andor_events, 'algebraic_net_time.xml')

        self.assertEqual(os.path.exists('algebraic_net_time.xml'), True)

        import libsbml
        r = libsbml.SBMLReader()
        d = r.readSBMLFromString('algebraic_net_time.xml')


    def test_stoichiometryPreservation(self):
        """ Test that products with stoichiometryMath set are output correctly. """
        SBMLInterface.toSBMLFile(algebraic_net, 'algebraic_net.xml')
        local_alg_net = SBMLInterface.fromSBMLFile('algebraic_net.xml')

        out_rxn = local_alg_net.reactions.get('out')
        self.assertEqual(out_rxn.product_stoichiometry['X1'], ['1 / 1'])

    def test_changeStoichiometry_SloppyCellModel(self):
        """ Test that products with stoichiometryMath set are output correctly. """
        SBMLInterface.toSBMLFile(algebraic_net, 'algebraic_net.xml')

        local_alg_net = copy.deepcopy(algebraic_net)
        out_rxn = local_alg_net.reactions.get('out')

        out_rxn.change_stoichiometry('T', -2)
        out_rxn.change_stoichiometry('X1', '2 / 2')
        
        self.assertEqual(out_rxn.reactant_stoichiometry, None)
        self.assertEqual(out_rxn.product_stoichiometry, None)
        self.assertEqual(out_rxn.stoichiometry['T'], -2)
        self.assertEqual(out_rxn.stoichiometry['X1'], '2 / 2')        

    def test_changeStoichiometry_SBMLModel(self):
        """ Test that products with stoichiometryMath set are output correctly. """
        SBMLInterface.toSBMLFile(algebraic_net, 'algebraic_net.xml')
        local_alg_net = SBMLInterface.fromSBMLFile('algebraic_net.xml')

        out_rxn = local_alg_net.reactions.get('out')

        out_rxn.change_stoichiometry('T', -2)
        out_rxn.change_stoichiometry('X1', '2 / 2')
        
        self.assertEqual(out_rxn.product_stoichiometry['T'], [-2])
        self.assertEqual(out_rxn.product_stoichiometry['X1'], ['2 / 2'])
        self.assertEqual(out_rxn.stoichiometry['T'], -2)
        self.assertEqual(out_rxn.stoichiometry['X1'], '2 / 2')

    def test_andor_funcs_toSBML(self):
        """ Test that events with complicated and_func and or_funcs are
        output to SBML correctly. """

        outfile = 'algebraic_net_andor_events.xml'

        piecewise_X0 = 'piecewise(0, or_func(gt(X0, 100), gt(X1, 100), lt(T, 0)), 1)'
        piecewise_X1 = 'piecewise(1, or_func(gt(X0, 100), gt(X1, 100), lt(T, 0)), 0)'
        logical_trigger = 'and_func(gt(1, 0), lt(1, 0), eq(10, 11))'

        

        if os.path.exists(outfile) == True:
            os.remove(outfile)        
        SBMLInterface.toSBMLFile(algebraic_net_andor_events, outfile)

        new_net = SBMLInterface.fromSBMLFile(outfile)
        logical_event = new_net.events.get('logical_event')

        self.assertEqual(piecewise_X0, logical_event.event_assignments.get('X0'))
        self.assertEqual(piecewise_X1, logical_event.event_assignments.get('X1'))
        self.assertEqual(logical_trigger, logical_event.trigger)

################################################################################
        
if _HAVE_SBML:
    suite = unittest.makeSuite(test_SBMLInterface)

if __name__ == '__main__':
    if _HAVE_SBML:
        unittest.main()
