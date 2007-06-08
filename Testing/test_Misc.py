import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_assignment
base_net = algebraic_net_assignment.copy()

class test_Misc(unittest.TestCase):
    def test_AssignedVarBug(self):
        """ Test handling of assigned variables initialized to concentration 
        'None'"""
        # This used to raise an exception.
        net = base_net.copy('test')
        net.add_species('tester', 'cell', None)
        net.add_assignment_rule('tester', 'X0')
        net.updateAssignedVars(1.0)

    def test_ChangingFunctionDefs(self):
        """
        Test whether changing function definitions are handled correctly.
        """
        net = Network('test')
        net.add_parameter('x', 0.0)
        net.add_rate_rule('x', 'f(1)')
        net.add_func_def('f', ['x'], 'x+2')
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.get_var_val('x', 10), 30)
        # It's not clear to me why this version wasn't causing failures
        # before...
        #net.remove_component('f')
        #net.add_func_def('f', ['x'], 'x+4')
        net.functionDefinitions.get('f').math = 'x+4'
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.get_var_val('x', 10), 50)

suite = unittest.makeSuite(test_Misc)
if __name__ == '__main__':
    unittest.main()
