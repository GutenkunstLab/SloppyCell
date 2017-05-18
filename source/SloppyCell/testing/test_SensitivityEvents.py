import unittest
import os

from SloppyCell.ReactionNetworks import *

import AlgTestNets
base_net = AlgTestNets.algebraic_net_assignment.copy()
base_net.remove_component('event0')
base_net.set_var_constant('k1', False)

dyn_var = 'X1'
opt_var = 'Keq'
eps = 1e-4
times = [0, 10]

def do_sens_and_fd(net):
    params = net.GetParameters()
    plus_params = params.copy()
    plus_params.set(opt_var, plus_params.get(opt_var) + eps)

    sens_traj = Dynamics.integrate_sensitivity(net, times, params, 
                                               fill_traj=False)

    traj_central = Dynamics.integrate(net, times, rtol=1e-10, params=params,
                                      fill_traj=False)
    traj_plus = Dynamics.integrate(net, times, rtol=1e-10, 
                                   params=plus_params, fill_traj=False)

    return sens_traj, traj_central, traj_plus

class test_SensitivityEvents(unittest.TestCase):
    def assert_over_all_vars(self, sens_traj, traj_central, traj_plus):
        for dyn_var in sens_traj.dynamicVarKeys:
            y_plus = traj_plus.get_var_traj(dyn_var)
            y_central = traj_central.get_var_traj(dyn_var)
            diff_sens = (y_plus-y_central)/eps

            sens_val = sens_traj.get_var_val_index((dyn_var, opt_var), -1)
            self.assertAlmostEqual(diff_sens[-1], sens_val, 4, 
                                   'failed for (%s, %s), %g vs %g'
                                   % (dyn_var, opt_var, diff_sens[-1],
                                      sens_val))

    def test_basic(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 4)', {'k1': 4.3})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_constant_assignment(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'lt(X0, Keq/4)', {'X1': 4.3})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_basic_delay(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2}, delay = 3)
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_chained_events(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'gt(time, 5)', 
                      {'k1': 2})
        net.add_event('event2', 'gt(k1, 1)', 
                      {'k1': 4})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

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
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_trigger_w_explicit_param(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'lt(X0, Keq/4)', {'k1': 4.3})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_assignment_w_explicit_param(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'lt(X0, Keq/4)', {'k1': 'X0/Keq'})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_logical_in_trigger_1(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'and_func(lt(X0, 0.4), gt(T, 0.3))', 
                      {'k1': '20+X0/Keq'})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_logical_in_trigger_2(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'and_func(lt(X0, 0.8), gt(T, 0.3))', 
                      {'k1': '20+X0/Keq'})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_piecewise_in_assignment_1(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'and_func(lt(X0, 0.8), gt(T, 0.3))', 
                      {'k1': 'piecewise(X0**Keq, lt(X0,T), X0/T + Keq)'})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)

    def test_piecewise_in_assignment_otherwise(self):
        net = base_net.copy('test')
        net.add_event('event 1', 'and_func(lt(X0, 0.8), gt(T, 0.3))', 
                      {'k1': 'piecewise(X0**Keq, gt(X0,T), X0/T + Keq)'})
        sens_traj, traj_central, traj_plus = do_sens_and_fd(net)
        self.assert_over_all_vars(sens_traj, traj_central, traj_plus)


suite = unittest.makeSuite(test_SensitivityEvents)
if __name__ == '__main__':
    unittest.main()
