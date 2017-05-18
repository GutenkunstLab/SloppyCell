import unittest

from SloppyCell.ReactionNetworks import *

import TestNetwork
base_net = TestNetwork.net.copy('test_logs')

class test_IntegrateWithLogs(unittest.TestCase):
    def test_basic(self):
        net = base_net.copy('test_basic')
        traj = Dynamics.integrate(net, [0, 10])
        net.integrateWithLogs = True
        log_traj = Dynamics.integrate(net, [0, 10])
        for var in net.variables.keys():
            norm_val = traj.get_var_val_index(var, -1)
            log_val = log_traj.get_var_val_index(var, -1)
            self.assertAlmostEqual(norm_val, log_val, 5, 'Failed for %s.' % var)

    def test_roots(self):
        net = base_net.copy('test_basic')
        net.add_event('test_event', 'lt(x, 0.5)', {'y': '3*x*y'})
        traj = Dynamics.integrate(net, [0, 10])
        net.integrateWithLogs = True
        log_traj = Dynamics.integrate(net, [0, 10])
        for var in net.variables.keys():
            norm_val = traj.get_var_val_index(var, -1)
            log_val = log_traj.get_var_val_index(var, -1)
            self.assertAlmostEqual(norm_val, log_val, 5, 'Failed for %s.' % var)

    def test_sens(self):
        net = base_net.copy('test_basic')
        net.add_event('test_event', 'lt(x, 0.5)', {'y': '3*x*y'})
        traj = Dynamics.integrate_sensitivity(net, [0, 10], rtol=1e-9)
        net.integrateWithLogs = True
        log_traj = Dynamics.integrate_sensitivity(net, [0, 10], rtol=1e-9)
        for dyn_var in net.dynamicVars.keys():
            for opt_var in net.optimizableVars.keys():
                norm_val = traj.get_var_val_index((dyn_var, opt_var), -1)
                log_val = log_traj.get_var_val_index((dyn_var, opt_var), -1)
                msg = 'Failed for (%s, %s) %f v %f'\
                        % (dyn_var, opt_var, norm_val, log_val)
                self.assertAlmostEqual(norm_val, log_val, 5, msg)
        
suite = unittest.makeSuite(test_IntegrateWithLogs)
if __name__ == '__main__':
    unittest.main()
