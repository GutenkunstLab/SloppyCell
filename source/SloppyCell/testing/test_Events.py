import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_assignment
base_net = algebraic_net_assignment.copy()
base_net.remove_component('event0')
base_net.set_var_constant('k1', False)

class test_Events(unittest.TestCase):
    def test_basic_delay(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2}, delay = 3)
        traj = Dynamics.integrate(net, [0, 5.5, 10])
        self.assertAlmostEqual(traj.get_var_val('k1', 5.5), 0.1)
        self.assertAlmostEqual(traj.get_var_val('k1', 10), 2)

    def test_chained_events(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2})
        net.add_event('event2', 'gt(k1, 1)', 
                      {'k1': 4})
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.get_var_val('k1', 10), 4)

    def test_longer_chained_events(self):
        net = base_net.copy('test')
        net.set_var_constant('k1', False)
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 1.5},
                      delay = 0)
        net.add_event('event2', 'gt(k1, 1)', 
                      {'k1': 2}, delay=0)
        net.add_event('event3', 'gt(k1, 1.75)', 
                      {'p_3': 1}, delay=1)
        traj = Dynamics.integrate(net, [0, 5.5, 6.5, 10])
        self.assertAlmostEqual(traj.get_var_val('k1', 10.0), 2)
        self.assertAlmostEqual(traj.get_var_val('p_3', 5.5), 0)
        self.assertAlmostEqual(traj.get_var_val('p_3', 6.5), 1)

    def test_piecewise_in_assignment(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(T, 0.3)', 
                      {'k1': 'piecewise(1, lt(X0,0.7), 2)'})
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.get_var_val('k1', 10), 1)

    def test_piecewise_in_assignment_otherwise(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(T, 0.3)', 
                      {'k1': 'piecewise(1, lt(X0,0.3), 2)'})
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.get_var_val('k1', 10), 2)

    def test_logical_in_trigger_1(self):
        net = base_net.copy('test')
        # This event should fire when T (dyn_var 2) gets to 0.3
        net.add_event('event 1', 'and_func(lt(X0, 0.8), gt(T, 0.3))', 
                      {'k1': 1})
        traj = Dynamics.integrate(net, [0, 10])

        self.assertAlmostEqual(traj.event_info[1][0][2], 0.3)

    def test_logical_in_trigger_2(self):
        net = base_net.copy('test')
        # This event should fire when X0 (dyn_var 2) gets to 0.8
        net.add_event('event 1', 'or_func(lt(X0, 0.8), gt(T, 0.3))', 
                      {'k1': 1})
        traj = Dynamics.integrate(net, [0, 10])
        self.assertAlmostEqual(traj.event_info[1][0][0], 0.8)

suite = unittest.makeSuite(test_Events)
if __name__ == '__main__':
    unittest.main()
