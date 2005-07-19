import copy
import unittest

import scipy

from TestNetwork import net, m, params

class test_ScaleFactors(unittest.TestCase):
    def test_fixed_scale_factors(self):
        """Test that fixed scale factors are handled correctly"""
        mCopy = copy.deepcopy(m)
        mCopy.get_expts()['expt1'].set_fixed_scale_factors({'x': 1.0,
                                            'y': 0.5})
        c = mCopy.cost(params)
        self.assertAlmostEqual(c, 139.3202824, 6,
                               'Failed on test of fixed scale factors.')
        
    def test_inconsistent_fixed_and_shared_scale_factors(self):
        """Test that inconsistent fixed and shared scale factors are flagged"""
        mCopy = copy.deepcopy(m)
        mCopy.get_expts()['expt1'].set_fixed_scale_factors({'x': 1.0,
                                            'y': 0.5})
        mCopy.get_expts()['expt1'].set_shared_scale_factors([['x', 'y']])
        self.assertRaises(ValueError, mCopy.cost, params)


    def test_fixed_and_shared_scale_factors(self):
        """Test that inconsistent fixed and shared scale factors are flagged"""
        mCopy = copy.deepcopy(m)
        mCopy.get_expts()['expt1'].set_shared_scale_factors([['x', 'y']])
        mCopy.get_expts()['expt1'].set_fixed_scale_factors({'x': 0.2})
        c = mCopy.cost(params)
        self.assertAlmostEqual(c, 182.61850579, 6,
                               'Failed on test of fixed shared scale factors.')

    def test_consistent_fixed_and_shared_scale_factors(self):
        """Test that consistent fixed and shared scale factors pass"""
        mCopy = copy.deepcopy(m)
        mCopy.get_expts()['expt1'].set_fixed_scale_factors({'x': 0.2,
                                            'y': 0.2})
        mCopy.get_expts()['expt1'].set_shared_scale_factors([['x', 'y']])
        c = mCopy.cost(params)
        self.assertAlmostEqual(c, 182.61850579, 6,
                               'Failed on test of dual fixed shared scale factors.')
    
    def test_shared_scale_factors(self):
        """Check that shared scale factors are computed correctly""" 
        mCopy = copy.deepcopy(m)
        mCopy.get_expts()['expt1'].set_shared_scale_factors([['x', 'y']])
        c = mCopy.cost(params)
        xsf = mCopy.internalVars['scaleFactors']['expt1']['x']
        self.assertAlmostEqual(xsf, 3.099782130, 6,
                                   'Shared scale factor computed incorrectly.')

suite = unittest.makeSuite(test_ScaleFactors, 'test')

if __name__ == '__main__':
    unittest.main()
