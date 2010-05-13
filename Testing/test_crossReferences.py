import unittest
import os
import copy

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *
import SloppyCell.daskr

from SloppyCell.daskr import daeint

# Load the fast reaction example from the SBML semantic test suite.
# To avoid extra dependencies on libsbml, we use verions built by SloppyCell.
from AlgTestNets import algebraic_net, algebraic_net_manual

tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 51)])

class test_crossReferences(unittest.TestCase):


    # Introducing two helper functions for the integration test, which is
    # repeated a few times below.
    def integration_tests(self, net):
        """ integrates the network to timepoints with known data which we test """
        traj = Dynamics.integrate(net, tlist_algebraic_net)
        self.data_test_cases(traj)

    def data_test_cases(self, traj):
        """ do the actual data comparison """
        self.assertAlmostEqual(traj.get_var_val('X0',4.8), 
                               0.618783392, 5)
        self.assertAlmostEqual(traj.get_var_val('X1',21.6), 
                               0.653837775, 5)
        self.assertAlmostEqual(traj.get_var_val('T', 29.6), 
                               0.138253942, 5)
        self.assertAlmostEqual(traj.get_var_val('S1', 40.0), 
                               0.018207409, 5)
        self.assertAlmostEqual(traj.get_var_val('S2', 16.8), 
                               0.210750878, 5)
        
    
    def test_manual_cross_refs(self):
        """ Test that a network with manual call to makeCrossReferences gets \
correct reults """

        net = algebraic_net_manual.copy()
        net._makeCrossReferences()

        self.integration_tests(net)


    def test_manual_without_makeCrossReferences(self):
        """ Test that after setting _manualCrossReferences_flag to True the \
compilation will fail if makeCrossReferences is NOT called """

        net = algebraic_net_manual.copy()
        self.assertRaises(AttributeError, Dynamics.integrate, net, tlist_algebraic_net)

    def test_manual_without_makeCrossReferences_2(self):
        """ Test that after setting _manualCrossReferences_flag to True the \
compilation will fail if makeCrossReferences is NOT called. Also test that after \
makeCrossReferences() is subsequently called the integration works."""

        net = algebraic_net_manual.copy()
        self.assertRaises(AttributeError, Dynamics.integrate, net, tlist_algebraic_net)
        net._makeCrossReferences()
        self.integration_tests(net)

    def test_manual_off(self):
        """ Test that after setting _manualCrossReferences_flag to True and then back \
integration still works """

        net = algebraic_net.copy()

        self.assertEqual(net._manualCrossReferences_flag, False)

        net._manualCrossReferences(flag=True)
        self.assertEqual(net._manualCrossReferences_flag, True)

        net._manualCrossReferences(flag=False)
        self.assertEqual(net._manualCrossReferences_flag, False)
        
        self.integration_tests(net)


    def test_manual_call_makeCrossReferences(self):
        """ Test that for a network where cross references are manual, once \
_makeCrossReferences is called then it will be automatic."""

        net = algebraic_net_manual.copy()

        self.assertEqual(algebraic_net_manual._manualCrossReferences_flag, True)
        self.assertEqual(net._manualCrossReferences_flag, True)

        net._makeCrossReferences()
        self.assertEqual(len(net.GetParameters()), 3)
        self.assertEqual(net._manualCrossReferences_flag, False)

        # Try adding a parameter to make sure the cross references are made automatically.
        net.addParameter('k_new', 1e-3)
        self.assertEqual(len(net.GetParameters()), 4)

        self.integration_tests(net)
        
        




################################################################################
        
suite = unittest.makeSuite(test_crossReferences)

if __name__ == '__main__':
    unittest.main()
