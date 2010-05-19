import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_assignment, algebraic_net
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

    def test_include_extra_info_in_event_info(self):
        """ Test that values of assigned variables can be included in event information """

        # First we test a trajectory that doesn't store the assigned
        # variables at events.
        net1 = base_net.copy('test')        
        net1.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2})
        net1.add_event('event2', 'gt(k1, 1)', 
                      {'k1': 4})
        traj1 = Dynamics.integrate(net1, [0, 5.5, 15])

        self.assertEqual(len(traj1.event_info), 3)
        
        # Then we test a trajectory that does store the assigned variables.
        net2 = net1.copy('test')        
        traj2 = Dynamics.integrate(net2, [0, 5.5, 15], include_extra_event_info =True)
        te, ye_pre, ye_post, ie, assign_e = traj2.event_info
        s_sum = assign_e[0]['S_sum']

        self.assertEqual(len(traj2.event_info), 5)
        self.assertAlmostEqual(s_sum, 0.29791465091325, 7)

    def test_trajectory_merge(self):
        """ Test basic trajectory merging """
        net = base_net.copy('test')
        times1 = scipy.linspace(0, 15, 100)
        times2 = scipy.linspace(15, 25, 100)        
        traj1 = Dynamics.integrate(net, times1, fill_traj=False)
        traj2 = Dynamics.integrate(net, times2, fill_traj=False)

        traj1.merge(traj2)

        times = traj1.get_times()
        self.assertEqual(len(times), 200)
        self.assertEqual(list(times), sorted(times))
        
    def test_merge_incompartible_trajs(self):
        """ Test the incompatible trajectories cause an exception to be raised """

        net1 = algebraic_net.copy()
        net2 = algebraic_net_assignment.copy()
        
        times1 = scipy.linspace(0, 15, 100)
        times2 = scipy.linspace(15, 25, 100)
        
        traj1 = Dynamics.integrate(net1, times1, fill_traj=False)
        traj2 = Dynamics.integrate(net2, times2, fill_traj=False)

        self.assertRaises(ValueError, traj1.merge, traj2)


        

suite = unittest.makeSuite(test_Trajectory)
if __name__ == '__main__':
    unittest.main()
