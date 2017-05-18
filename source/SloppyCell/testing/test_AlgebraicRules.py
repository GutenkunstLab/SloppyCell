import unittest
import os

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *
import SloppyCell.daskr

from SloppyCell.daskr import daeint

# Load the fast reaction example from the SBML semantic test suite.
# To avoid extra dependencies on libsbml, we use verions built by SloppyCell.
from AlgTestNets import algebraic_net, algebraic_net_assignment, \
                         algebraic_net_multi, algebraic_net_under
tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 51)])

class test_AlgebraicRules(unittest.TestCase):
    def test_basic(self):
        """ Basic test of Algebraic Rules """
        algebraic_traj = Dynamics.integrate(algebraic_net, tlist_algebraic_net)

        self.assertAlmostEqual(algebraic_traj.get_var_val('X0',4.8), 
                               0.618783392, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('X1',21.6), 
                               0.653837775, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('T', 29.6), 
                               0.138253942, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S1', 40.0), 
                               0.018207409, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S2', 16.8), 
                               0.210750878, 5)

    def test_several_algebraic_rules(self):
        """ Test of algebraic system with 3 algebraic rules """

        algebraic_traj = Dynamics.integrate(algebraic_net_multi, tlist_algebraic_net, redirect_msgs=False)

        self.assertAlmostEqual(algebraic_traj.get_var_val('X0',4.8), 
                               0.618783392, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('X1',21.6), 
                               0.653837775, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('T', 29.6), 
                               0.138253942, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S1', 40.0), 
                               0.018207409, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S2', 16.8), 
                               0.210750878, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S3', 16.8), 
                               0.1, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S4', 16.8), 
                               0.15, 5)

    def test_assignment_in_algebraic(self):
        """ Test that algebraic rhs functions that are functions of assignment
            rules are parsed correctly"""

        # integrate using the same time points as the base case
        algebraic_traj = Dynamics.integrate(algebraic_net_assignment, 
                                            tlist_algebraic_net, rtol=1e-9)

        # make sure that the correct variables were identified as algebraic

        alg_vars = algebraic_net_assignment.algebraicVars
    
        self.assertEqual(alg_vars.has_key('S1'), True)
        self.assertEqual(alg_vars.has_key('T'), False)
        self.assertEqual(alg_vars.has_key('S_sum'), False)
        self.assertEqual(alg_vars.has_key('Alg_Rule_RHS'), False)

        # p_1 is a paramter that was set to be non-constant, but which is
        # implicitly constant because no equations modify it.  We need to check
        # to make sure it's not being identified as an algebraic variable.
        self.assertEqual(alg_vars.has_key('p_1'), False)

        # p_2 is modified by an event assignment, so it should not be identified
        # as an algebraic variable.
        self.assertEqual(alg_vars.has_key('p_2'), False)
        

        # finally, we check to make sure the algebraic variables list is the
        # correct length
        self.assertEqual(len(alg_vars), 1)

        self.assertAlmostEqual(algebraic_traj.get_var_val('X0', 4.8),
                               0.618783392, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('X1', 21.6), 
                               0.653837775, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('T', 29.6), 
                               0.138253942, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S1', 40.0), 
                               0.018207409, 5)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S2', 16.8), 
                               0.210750878, 5)


    def test_low_tolerances(self):
        """ Test that algebraic rules still work correctly when low tolerances
        are enforced"""

        reltol = scipy.array([1e-10, 1e-10, 1e-10, 1e-10])

        algebraic_traj = Dynamics.integrate(algebraic_net, tlist_algebraic_net,
                                              reltol)

        self.assertAlmostEqual(algebraic_traj.get_var_val('X0', 4.8),
                               0.618783392, 8)
        self.assertAlmostEqual(algebraic_traj.get_var_val('X1', 21.6), 
                               0.653837775, 8)
        self.assertAlmostEqual(algebraic_traj.get_var_val('T', 29.6), 
                               0.138253942, 8)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S1', 40.0), 
                               0.018207409, 8)
        self.assertAlmostEqual(algebraic_traj.get_var_val('S2', 16.8),
                               0.210750878, 8)

    def test_log_integration(self):
        """ Test that algebraic rules still work correctly when integration
        is done in terms of log chem values."""

        reltol = scipy.array([1e-10, 1e-10, 1e-10, 1e-10])

        algebraic_traj = Dynamics.integrate(algebraic_net, tlist_algebraic_net,
                                            reltol, fill_traj=False)

        t0 = tlist_algebraic_net[1]
        vals = algebraic_traj.get_var_vals(t0)

        net = algebraic_net.copy('log_test')
        net.set_var_ics(vals)
        net.integrateWithLogs = True
        log_traj = Dynamics.integrate(net, tlist_algebraic_net,
                                      reltol, fill_traj=False)

        self.assertAlmostEqual(log_traj.get_var_val('X0', 4.8-t0),
                               0.618783392, 8)
        self.assertAlmostEqual(log_traj.get_var_val('X1', 21.6-t0), 
                               0.653837775, 8)
        self.assertAlmostEqual(log_traj.get_var_val('T', 29.6-t0), 
                               0.138253942, 8)
        self.assertAlmostEqual(log_traj.get_var_val('S1', 40.0-t0), 
                               0.018207409, 8)
        self.assertAlmostEqual(log_traj.get_var_val('S2', 16.8-t0),
                               0.210750878, 8)

    def test_underdetermined(self):
        """ Test that an underdetermined system raises a ValueError """
        
        # algebraic_net_under is underdetermined (no rule for determining parameter P)
        # check that system raises ValueError in compilation
        self.assertRaises((ValueError), Dynamics.integrate,
                          algebraic_net_under, tlist_algebraic_net)



################################################################################
        
suite = unittest.makeSuite(test_AlgebraicRules)

if __name__ == '__main__':
    unittest.main()
