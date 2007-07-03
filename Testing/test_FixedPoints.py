import unittest

import scipy

from SloppyCell.ReactionNetworks import *

lorenz = Network('lorenz')
lorenz.add_compartment('basic')
lorenz.add_species('x', 'basic', 0.5)
lorenz.add_species('y', 'basic', 0.5)
lorenz.add_species('z', 'basic', 0.5)
lorenz.add_parameter('sigma', 1.0)
lorenz.add_parameter('r', 2.0)
lorenz.add_parameter('b', 2.0)
lorenz.add_rate_rule('x', 'sigma*(y-x)')
lorenz.add_rate_rule('y', 'r*x - y - x*z')
lorenz.add_rate_rule('z', 'x*y - b*z')

class test_fixedpoints(unittest.TestCase):
    def test_basic(self):
        """ Test basic fixed-point finding """
        net = lorenz.copy('test')
        fp = Dynamics.dyn_var_fixed_point(net, dv0=[1,1,1], with_logs=False)
        # This should find the fixed-point [sqrt(2), sqrt(2), 1]
        self.assertAlmostEqual(fp[0], scipy.sqrt(2), 6, 'Failed on basic 1,0.')
        self.assertAlmostEqual(fp[1], scipy.sqrt(2), 6, 'Failed on basic 1,1.')
        self.assertAlmostEqual(fp[2], 1, 6, 'Failed on basic 1,2.')

        fp = Dynamics.dyn_var_fixed_point(net, dv0=[-0.1,-0.1,-0.1], 
                                          with_logs=False)
        # This should find the fixed-point [0, 0, 0]
        self.assertAlmostEqual(fp[0], 0, 6, 'Failed on basic 2,0.')
        self.assertAlmostEqual(fp[1], 0, 6, 'Failed on basic 2,1.')
        self.assertAlmostEqual(fp[2], 0, 6, 'Failed on basic 2,2.')

    def test_withlogs(self):
        """ Test fixed-point finding with logs """
        net = lorenz.copy('test')
        fp = Dynamics.dyn_var_fixed_point(net, dv0=[1,1,1], with_logs=True)
        # This should find the fixed-point [sqrt(2), sqrt(2), 1]
        self.assertAlmostEqual(fp[0], scipy.sqrt(2), 6, 'Failed on logs 1,0.')
        self.assertAlmostEqual(fp[1], scipy.sqrt(2), 6, 'Failed on logs 1,1.')
        self.assertAlmostEqual(fp[2], 1, 6, 'Failed on logs 1,2.')

        fp = Dynamics.dyn_var_fixed_point(net, dv0=[0.1,0.1,0.1], 
                                          with_logs=True)
        # This should find the fixed-point [0, 0, 0]
        self.assertAlmostEqual(fp[0], 0, 6, 'Failed on logs 2,0.')
        self.assertAlmostEqual(fp[1], 0, 6, 'Failed on logs 2,1.')
        self.assertAlmostEqual(fp[2], 0, 6, 'Failed on logs 2,2.')

    def test_stability(self):
        net = lorenz.copy('test')
        # The sqrt(b*(r-1)), sqrt(b*(r-1)), r-1 fixed point is stable for r < rH
        #  Strogatz, Nonlinear Dynamics and Chaos (p. 316)
        fp, stable = Dynamics.dyn_var_fixed_point(net, dv0=[1,1,1], 
                                                  stability=True)
        self.assertEqual(stable, -1, 'Failed to classify stable fixed point')

        # (0,0,0) is a saddle here
        fp, stable = Dynamics.dyn_var_fixed_point(net, dv0=[0.01,0.01,0.01], 
                                                  stability=True)
        self.assertEqual(stable, 0, 'Failed to classify saddle')

        # (0,0,0) is a stable node here
        net.set_var_ic('r', 0.5)
        fp, stable = Dynamics.dyn_var_fixed_point(net, dv0=[0.1,0.1,0.1], 
                                                  stability=True)
        self.assertEqual(stable, -1, 'Failed to classify stable fixed point')

        # Now make the far fixed point a saddle...
        net.set_var_ic('sigma', 6.0)
        net.set_var_ic('r', 25)
        fp, stable = Dynamics.dyn_var_fixed_point(net, dv0=[10,10,10], 
                                                  stability=True)
        self.assertEqual(stable, 0, 'Failed to classify saddle')

suite = unittest.makeSuite(test_fixedpoints)

if __name__ == '__main__':
    unittest.main()
