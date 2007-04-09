import unittest
import os

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *


# Load a similar algebraic rule where assignment variables are used in the rule.
file_path = os.path.join('SBML_files',
                         'algebraicRules-assignment_in_algebraic.xml')
base_net = IO.from_SBML_file(file_path,
                             'algebraicRules-assignment_in_algebraic')


class test_Interpolation(unittest.TestCase):
    def test_basic(self):
        """ Test that interpolated trajectories don't crash. """
        traj = Dynamics.integrate(base_net, [0, 10], fill_traj=True)
        traj.build_interpolated_traj()
    
        traj.evaluate_interpolated_traj('X1', [1, 2])
        traj.evaluate_interpolated_traj('X1', 1)
        traj.evaluate_interpolated_traj('X1', [1, 2], der=3)
        
suite = unittest.makeSuite(test_Interpolation)

if __name__ == '__main__':
    unittest.main()
