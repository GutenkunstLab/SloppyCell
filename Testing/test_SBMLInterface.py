import unittest
import os

import scipy
import SloppyCell.Utility as Utility
from SloppyCell.ReactionNetworks import *

# Load the fast reaction example from the SBML semantic test suite.
# To avoid extra dependencies on libsbml, we use verions built by SloppyCell.
from AlgTestNets import algebraic_net, algebraic_net_assignment, algebraic_net_multi
tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 51)])

class test_SBMLInterface(unittest.TestCase):
    def test_tofile(self):
        """ Basic test of SBML output """

        # First remove the output file in case it was produced in a
        # previous test
        if os.path.exists('algebraic_net.xml') == True:
            os.remove('algebraic_net.xml')

        # Now output the SBML and make sure it got created        
        SBMLInterface.toSBMLFile(algebraic_net, 'algebraic_net.xml')

        self.assertEqual(os.path.exists('algebraic_net.xml'), True)


################################################################################
        
suite = unittest.makeSuite(test_SBMLInterface)

if __name__ == '__main__':
    unittest.main()
