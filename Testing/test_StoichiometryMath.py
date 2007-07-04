import unittest
import os

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *


# Load the stoichiometryMathML example from the SBML semantic test suite.
# To avoid extra dependencies on libsbml, we use verions built by SloppyCell.
from StoichTestNets import stoichMath_net

tlist_stoichMath_net = scipy.array([0] + [0.04*x for x in range(1, 51)])

class test_StoichiometryMath(unittest.TestCase):
    def test_basic(self):
        """ Basic test of stoichiometry math """
        stoichMath_traj = Dynamics.integrate(stoichMath_net, tlist_stoichMath_net)

        self.assertAlmostEqual(stoichMath_traj.get_var_val('A',0.32), 
                               0.62283185811441, 5)
        self.assertAlmostEqual(stoichMath_traj.get_var_val('A',1.28), 
                               0.0129772159879822, 5)
        self.assertAlmostEqual(stoichMath_traj.get_var_val('A', 1.96), 
                               0.000623972556903647, 5)
        self.assertAlmostEqual(stoichMath_traj.get_var_val('B', 0.16), 
                               0.153505642887153, 5)
        self.assertAlmostEqual(stoichMath_traj.get_var_val('B', 1.84), 
                               0.44697495892817, 5)


################################################################################
        
suite = unittest.makeSuite(test_StoichiometryMath)

if __name__ == '__main__':
    unittest.main()
