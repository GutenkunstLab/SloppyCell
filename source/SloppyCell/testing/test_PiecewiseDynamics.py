import unittest

from SloppyCell.ReactionNetworks import *

import TestNetwork

class test_PiecewiseDynamics(unittest.TestCase):
    #XXX: Assignment rules currently not supported. To do so, add a vector
    #     version of piecewise to Trajectory_mod.py
    #def test_assignment(self):
    #    net = TestNetwork.net.copy('piecewise')
    #    net.disable_deriv_funcs()
    #    net.compile(disable_c=True)

    #    net.addSpecies('C', 'basic')
    #    net.add_assignment_rule('C', 'piecewise(2.0, time < 1.0, 1.0)')
    #    traj = Dynamics.integrate(net, [0, 0.5, 1.5])

    def test_raterule(self):
        net = TestNetwork.net.copy('piecewise')
        net.addSpecies('C', 'basic', 0)
        net.add_rate_rule('C', 'piecewise(2.0, time < 1.0, 1.0)')

        net.disable_deriv_funcs()
        net.disable_c = True

        traj = Dynamics.integrate(net, [0, 0.5, 1.5])
        self.assertAlmostEqual(traj.get_var_val('C', 1.5), 2.5, 3)

suite = unittest.makeSuite(test_PiecewiseDynamics)
if __name__ == '__main__':
    unittest.main()
