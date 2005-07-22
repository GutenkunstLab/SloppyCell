import unittest

import os,copy
import SloppyCell.ReactionNetworks.SBMLInterface as SBMLInt
from TestNetwork import net
net = copy.deepcopy(net)
net.compile()
 
class test_SBMLFunctions(unittest.TestCase):
    #def setUp(self) :
#	self.net = net    
    
#    def tearDown(self) :
#    	self.net = None
 
    def test_ToAndFromSBMLFile(self):
        """Testing the reading and writing of an SBML file from and to a network"""
	SBMLInt.toSBMLFile(net,'netSBML')
	netread = SBMLInt.fromSBMLFile('netSBML')
	os.system('rm netSBML')	

suite = unittest.makeSuite(test_SBMLFunctions, 'test')

if __name__ == '__main__':
    unittest.main()
