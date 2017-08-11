"""
Tests for parallel code
"""
import copy, scipy, unittest

from SloppyCell.ReactionNetworks import *

import TestNetwork
m2 = copy.deepcopy(TestNetwork.m2)
params = m2.get_params()

class test_parallel(unittest.TestCase):
    def test_cost(self):
        """ Test basic cost evaluation """
        self.assertAlmostEqual(m2.cost(params), 57.629, 2)

    # It's difficult to test exception handling, because I want to change
    #  the networks, but our synchronization is really primative.

    #def test_cost_exceptions_first(self):
    #    """ Test cost exceptions in first network"""
    #    m3 = copy.deepcopy(m2)
    #    temp = m3.get_calcs()['test'].get_ddv_dt
    #    # This should raise an exception
    #    m3.get_calcs()['test'].get_ddv_dt = lambda x,t: 'a'
    #    self.assertRaises(Exception, m3.cost, params)
    #    # But the next cost evaluation should still be correct
    #    m3.get_calcs()['test'].get_ddv_dt = temp
    #    self.assertAlmostEqual(m3.cost(params), 57.629, 2)

    #def test_cost_exceptions_second(self):
    #    """ Test cost exceptions in second network"""
    #    m3 = copy.deepcopy(m2)
    #    temp = m3.get_calcs()['test'].get_ddv_dt
    #    # This should raise an exception that gets passed to the master node
    #    m3.get_calcs()['test2'].get_ddv_dt = lambda x,t: 'a'
    #    self.assertRaises(Exception, m3.cost, params)
    #    # But the next cost evaluation should still be correct
    #    m3.get_calcs()['test2'].get_ddv_dt = temp
    #    self.assertAlmostEqual(m3.cost(params), 57.629, 2)

    def test_JtJ(self):
        """ Test that JtJ calculation doesn't crash """
        jtj = m2.GetJandJtJInLogParameters(scipy.log(params))

if num_procs > 1:
    suite = unittest.makeSuite(test_parallel)

if __name__ == '__main__':
    if num_procs == 1:
        print 'Only one processor detected! Not running in parallel!'
