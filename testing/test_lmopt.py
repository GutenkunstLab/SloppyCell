import copy
import unittest

import scipy

from SloppyCell.ReactionNetworks import *
import SloppyCell
import TestNetwork
import SloppyCell.lmopt as lmopt
import os

net = TestNetwork.net
m = copy.deepcopy(TestNetwork.m)
p = m.params.__copy__()

class test_lmopt(unittest.TestCase):
    def test_lm_fmin(self):
        """Test that Levenberg-Marquardt algortihm runs correctly"""
        p = m.params.__copy__()
        p.update([1.0,1.0])
        res_func = m.res_log_params
        def res_func_prime(x) :
            j,jtj = m.GetJandJtJInLogParameters(x)
            jarray = scipy.zeros((len(j.keys()),len(m.params)),scipy.float_)
            for resind, resname in enumerate(m.residuals.keys()) :
                jarray[resind,:] = j.get(resname)
            return jarray
        
        #print "\n Initial cost", m.cost(m.params)
        pbest = lmopt.fmin_lm(res_func,
                scipy.log(p),res_func_prime,maxiter=20,disp=0,
                full_output=1) 
        m.params.update(scipy.exp(pbest[0]))
        newcost = m.cost(m.params) 
        #print "Cost after 20 iterations", newcost
        # The intial cost is about 25, the final cost should be 
        # about 1.0e-6
        self.assertAlmostEqual(newcost, 0.0, 4, 'Failed on lm minimum')
    
    def test_lm_fminNoJ(self):
        """Test that Levenberg-Marquardt algortihm which just uses JtJ 
        runs correctly"""
        p = m.params.__copy__()
        p.update([1.0,1.0])
        cost_func = m.cost_log_params
        def grad_and_lmhess_func(x) :
            j,jtj = m.GetJandJtJInLogParameters(x)
            jarray = scipy.zeros((len(j.keys()),len(m.params)),scipy.float_)
            for resind, resname in enumerate(m.residuals.keys()) :
                jarray[resind,:] = j.get(resname)
            resvals = m.res_log_params(x)
            grad = scipy.dot(scipy.transpose(jarray),resvals)
            return grad,jtj
        
        #print "\n Initial cost", m.cost(p)
        pbest = lmopt.fmin_lmNoJ(cost_func,
                scipy.log(p),grad_and_lmhess_func,maxiter=20,disp=0,
                full_output=1) 
        m.params.update(scipy.exp(pbest[0]))
        newcost = m.cost(m.params) 
        #print "Cost after 20 iterations", newcost
        # The intial cost is about 25, the final cost should be 
        # about 1.0e-6
        self.assertAlmostEqual(newcost, 0.0, 4, 'Failed on lm NoJ minimum')

    def test_lm_fmin_scale(self):
        """Test that the scale invariant Levenberg-Marquardt algortihm 
        runs correctly"""
        p = m.params.__copy__()
        p.update([1.0,1.0])
        res_func = m.res_log_params
        def res_func_prime(x) :
            j,jtj = m.GetJandJtJInLogParameters(x)
            jarray = scipy.zeros((len(j.keys()),len(m.params)),scipy.float_)
            for resind, resname in enumerate(m.residuals.keys()) :
                jarray[resind,:] = j.get(resname)
            return jarray
        
        #print "\n Initial cost", m.cost(m.params)
        pbest = lmopt.fmin_lm_scale(res_func,
                scipy.log(p),res_func_prime,maxiter=20,disp=0,
                full_output=1) 
        m.params.update(scipy.exp(pbest[0]))
        newcost = m.cost(m.params) 
        #print "Cost after 20 iterations", newcost
        # The intial cost is about 25, the final cost should be 
        # about 1.3 (this algorithm is not as efficient as the 
        # previous two, it's use is mainly for problems where different
        # parameters are on completely different scales)
        #self.assertLessThan(newcost, 2.0,'Failed on lm scale invariant minimum')

suite = unittest.makeSuite(test_lmopt)

if __name__ == '__main__':
    unittest.main()
