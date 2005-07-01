import unittest

import scipy

from SloppyCell import *

from TestNetwork import net

class test_ScaleFactors(unittest.TestCase):
    def setUp(self):
        self.params = [1.0, 2.0]
    
        # Make the experiment
        self.expt1 = Experiment('expt1')
        data = {'x': {1.0: (1.0, 0.1),
                      2.0: (1.0, 0.1)
                      },
                'y': {0.5: (2.0, 0.2),
                      2.5: (0.5, 0.05)
                      }
                     }
        self.expt1.set_data({'test': data})

        self.m = Model(ExperimentCollection([self.expt1]), 
                       CalculationCollection([net]))
        
    def test_fixed_scale_factors(self):
        """Test that fixed scale factors are handled correctly"""
        self.expt1.set_fixed_scale_factors({'x': 1.0,
                                            'y': 0.5})
        c = self.m.cost(self.params)
        self.assertAlmostEqual(c, 278.64056486, 6,
                                   'Failed on test of fixed scale factors.')
        
    def test_inconsistent_fixed_and_shared_scale_factors(self):
        """Test that inconsistent fixed and shared scale factors are flagged"""
        self.expt1.set_fixed_scale_factors({'x': 1.0,
                                            'y': 0.5})
        self.expt1.set_shared_scale_factors([['x', 'y']])
        self.assertRaises(ValueError, self.m.cost, self.params)


    def test_fixed_and_shared_scale_factors(self):
        """Test that inconsistent fixed and shared scale factors are flagged"""
        self.expt1.set_shared_scale_factors([['x', 'y']])
        self.expt1.set_fixed_scale_factors({'x': 0.2})
        c = self.m.cost(self.params)
        self.assertAlmostEqual(c, 365.237011585, 6,
                                   'Failed on test of fixed shared scale factors.')

    def test_consistent_fixed_and_shared_scale_factors(self):
        """Test that consistent fixed and shared scale factors pass"""
        self.expt1.set_fixed_scale_factors({'x': 0.2,
                                            'y': 0.2})
        self.expt1.set_shared_scale_factors([['x', 'y']])
        c = self.m.cost(self.params)
        self.assertAlmostEqual(c, 365.237011585, 6,
                                   'Failed on test of dual fixed shared scale factors.')
    
    def test_shared_scale_factors(self):
        """Check that shared scale factors are computed correctly""" 
        self.expt1.set_shared_scale_factors([['x', 'y']])
        c = self.m.cost(self.params)
        xsf = self.m.internalVars['scaleFactors']['expt1']['x']
        self.assertAlmostEqual(xsf, 3.099782130, 6,
                                   'Shared scale factor computed incorrectly.')

suite = unittest.makeSuite(test_ScaleFactors, 'test')

if __name__ == '__main__':
    unittest.main()
