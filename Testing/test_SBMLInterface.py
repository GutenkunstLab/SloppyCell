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
from AlgTestNets import algebraic_net, algebraic_net_assignment, algebraic_net_multi
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


################################################################################
        
if _HAVE_SBML:
    suite = unittest.makeSuite(test_SBMLInterface)

if __name__ == '__main__':
    if _HAVE_SBML:
        unittest.main()
