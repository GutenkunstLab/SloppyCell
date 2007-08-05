import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_assignment
base_net = algebraic_net_assignment.copy()
base_net.remove_component('event0')
base_net.set_var_constant('k1', False)

class test_Trajectory(unittest.TestCase):
    def test_interpolate_basic(self):
        net = base_net.copy('test')        
        traj = Dynamics.integrate(net, [0, 5.5, 15])
        traj.build_interpolated_traj()
        tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 15)])
        times, interp_traj =  traj.evaluate_interpolated_traj('X0', tlist_algebraic_net)

        self.assertAlmostEqual(interp_traj[1], 0.923116346390738, 5)
        self.assertAlmostEqual(interp_traj[4], .726149037073939, 5)
        self.assertAlmostEqual(interp_traj[6], 0.618783391805766, 5)
        self.assertAlmostEqual(interp_traj[8], 0.527292424043448, 5)
        self.assertAlmostEqual(interp_traj[10], 0.449328964119144, 5)
        self.assertAlmostEqual(interp_traj[14], 0.326279794615302, 5)

    def test_interpolate_chained_event(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2})
        net.add_event('event2', 'gt(k1, 1)', 
                      {'k1': 4})
        traj = Dynamics.integrate(net, [0, 5.5, 10, 15])
        self.assertAlmostEqual(traj.get_var_val('k1', 10), 4.0)
        
        traj.build_interpolated_traj()
        tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 15)])
        times, interp_traj =  traj.evaluate_interpolated_traj('X0', tlist_algebraic_net)

        self.assertAlmostEqual(interp_traj[1], 0.923116346390738, 5)
        self.assertAlmostEqual(interp_traj[4], .726149037073939, 5)
        self.assertAlmostEqual(interp_traj[6], 0.618783391805766, 5)
        self.assertAlmostEqual(interp_traj[8], 2.24340106e-03, 5)
        self.assertAlmostEqual(interp_traj[10], 4.37991652e-06, 5)
        self.assertAlmostEqual(interp_traj[14], 3.76502117e-12, 5)
        

suite = unittest.makeSuite(test_Trajectory)
if __name__ == '__main__':
    unittest.main()
