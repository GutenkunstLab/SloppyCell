import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_assignment
base_net = algebraic_net_assignment.copy()

class test_Misc(unittest.TestCase):
    def test_AssignedVarBug(self):
        """ Test handling of assigned variables initialized to concentration 
        'None'"""
        net = base_net.copy('test')
        net.add_species('tester', 'cell', None)
        net.add_assignment_rule('tester', 'X0')
        net.updateAssignedVars(1.0)

suite = unittest.makeSuite(test_Misc)
if __name__ == '__main__':
    unittest.main()
