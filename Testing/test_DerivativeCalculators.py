import copy
import unittest

import scipy

from SloppyCell.ReactionNetworks import *

import TestNetwork

net = TestNetwork.net
m = TestNetwork.m

class test_DerivativeCalculators(unittest.TestCase):
    def test_GetJandJtJ(self):
        """Test that finite difference Jacobian matches sensitivity Jacobian"""
       	p = m.params.__copy__()
	j,jtj = m.GetJandJtJ(p)

	# compare with finite diffences
	net.fdSensitivities = True
	p = m.params.__copy__()	
	jfd,jtjfd = m.GetJandJtJ(p) 
       	maxdiff = 0.0	
	for resname in j.keys() :
		avg = .5*abs(scipy.asarray(j[resname]))+.5*abs(scipy.asarray(jfd[resname]))	
		diff = abs(scipy.asarray(j[resname])-scipy.asarray(jfd[resname]))
		for i in range(len(diff)) :
			if avg[i] > 1.0 :  # use relative error for big entries	
				diff[i] = diff[i]/avg[i]
		diff = max(diff)	
		if diff > maxdiff :
			maxdiff = diff
	self.assertAlmostEqual(maxdiff, 0.0, 3,'Failed on J accuracy')
	
	maxdiff = 0.0	
	for i in range(len(p)) :
		for j in range(len(p)) :
			avg = .5*abs(jtj[i][j])+abs(jtjfd[i][j])
			diff = abs(jtj[i][j]-jtjfd[i][j])
			if avg > 1.0 : 	
				diff = diff/avg	
			if diff > maxdiff :
				maxdiff = diff	 
	self.assertAlmostEqual(maxdiff, 0.0, 3, 'Failed on JtJ accuracy')

    def test_HessianFunctions(self) :
	""" Test CalcHessian against CalcHessianUsingResiduals with and without log parameters """
	p = m.params.__copy__()
	h = m.CalcHessianUsingResiduals(p,1.0e-6,moreAcc=False)
	p = m.params.__copy__()	
	hfd = m.CalcHessian(p,1.0e-5)
	maxdiff = 0.0

	for i in range(len(p)) :
                for j in range(i,len(p)) :
                       	avg = .5*abs(h[i][j]) + .5*abs(hfd[i][j]) 
			diff = abs(h[i][j]-hfd[i][j])
			if avg > 1.0 :
				diff = diff/avg 
			if diff > maxdiff :
                                maxdiff = diff
        self.assertAlmostEqual(maxdiff, 0.0, 1, 'Failed on Hessian accuracy') # within 10% error is fine
	p = m.params.__copy__()
	h = m.CalcHessianUsingResidualsInLogParams(scipy.log(p),1.0e-6,moreAcc=False)
	p = m.params.__copy__()	
	hfd = m.CalcHessianInLogParameters(p,1.0e-5)
	maxdiff = 0.0
        for i in range(len(p)) :
                for j in range(len(p)) :
                       	avg = .5*abs(h[i][j]) + .5*abs(hfd[i][j]) 
			diff = abs(h[i][j]-hfd[i][j])
                       	if avg > 1.0 :
                                diff = diff/avg 
			if diff > maxdiff :
                                maxdiff = diff
        self.assertAlmostEqual(maxdiff, 0.0, 1, 'Failed on Hessian in log parameters accuracy') # within 10% error is fine
	 
suite = unittest.makeSuite(test_DerivativeCalculators)

if __name__ == '__main__':
    unittest.main()
