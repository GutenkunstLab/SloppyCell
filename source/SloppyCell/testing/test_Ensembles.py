import unittest

import scipy

from SloppyCell.ReactionNetworks import *

class test_Ensembles(unittest.TestCase):
    def test_PCA(self):
        """
        Test Principle Components Analysis function
        """
        # Generate gaussian set of log parameters with some variation
        #  in how tightly constrained it is
        log_p = scipy.randn(1000, 10)
        log_p[:,4] /= 2.0
        log_p[:,8] *= 5.0

        u,v = Ensembles.PCA_eig_log_params(scipy.exp(log_p))
        # Stiffest direction
        self.assertAlmostEqual(u[0], 4, 0)         # Eigenvalue
        self.assertAlmostEqual(abs(v[4,0]), 1, 2)  # Eigenvector

        # Sloppiest direction
        self.assertAlmostEqual(u[-1], 1/25., 2)
        self.assertAlmostEqual(abs(v[8,-1]), 1, 2)

        # One in the middle
        self.assertAlmostEqual(u[3], 1., 0)


suite = unittest.makeSuite(test_Ensembles)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
