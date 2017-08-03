import copy
import unittest

import scipy

from TestNetwork import net
net = copy.deepcopy(net)
net.setInitialVariableValue('A', 1.0)
net.setInitialVariableValue('B', 2.0)

class test_ics(unittest.TestCase):
    def test_default_initial_conditions(self):
        """Test that default ICs are handled correctly"""
        test_net = net.copy()
        traj = test_net.integrate(scipy.linspace(0, 5, 5))
        ICx = traj.getVariableTrajectory('x')[0]
        ICy = traj.getVariableTrajectory('y')[0]
        self.assertAlmostEqual(ICx, 1.0, 6, 'Failed on default IC')
        self.assertAlmostEqual(ICy, 2.0, 6, 'Failed on default IC')

    def test_resetting_initial_conditions(self):
        """Test resetting of ICs"""
        test_net = net.copy()
        test_net.set_initial_var_value('x', 0.5)
        traj = test_net.integrate(scipy.linspace(0, 5, 5))
        ICx = traj.getVariableTrajectory('x')[0]
        self.assertAlmostEqual(ICx, 0.5, 6, 'Failed on resetting IC')

    def test_parameter_ics(self):
        """Test parameters as ICs"""
        test_net = net.copy()
        test_net.set_initial_var_value('x', 'A')
        traj = test_net.integrate(scipy.linspace(0, 5, 5))
        ICx = traj.getVariableTrajectory('x')[0]
        self.assertAlmostEqual(ICx, 1.0, 6, 'Failed on parameter IC')

    def test_resetting_parameter(self):
        """Test changing parameters as ICs"""
        test_net = net.copy()
        test_net.set_initial_var_value('x', 'A')
        test_net.set_initial_var_value('A', 0.9)
        traj = test_net.integrate(scipy.linspace(0, 5, 5))
        ICx = traj.getVariableTrajectory('x')[0]
        self.assertAlmostEqual(ICx, 0.9, 6, 'Failed on changing parameter IC')

    def test_expression_ICs(self):
        """Test math expression as IC"""
        test_net = net.copy()
        test_net.set_initial_var_value('x', 'A + 1.5*B')
        traj = test_net.integrate(scipy.linspace(0, 5, 5))
        ICx = traj.getVariableTrajectory('x')[0]
        self.assertAlmostEqual(ICx, 4.0, 6, 'Failed on changing parameter IC')


suite = unittest.makeSuite(test_ics)

if __name__ == '__main__':
    unittest.main()
