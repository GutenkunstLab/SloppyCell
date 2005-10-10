import copy
import os
import unittest

try:
    import SloppyCell.ReactionNetworks.SBMLInterface as SBMLInt
    _HAVE_SBML = True
except ImportError:
    _HAVE_SBML = False

from TestNetwork import net
net = copy.deepcopy(net)
net.compile()
 
class test_SBMLFunctions(unittest.TestCase):
    def test_ToAndFromSBMLFile(self):
        """Testing the reading and writing of an SBML file from and to a network"""
        filename = 'net.sbml'
	SBMLInt.toSBMLFile(net, filename)
	netread = SBMLInt.fromSBMLFile(filename)
	os.remove(filename)	

no_libsbml_msg = 'libsbml not installed. SBML import and export will not be available!'
if _HAVE_SBML:
    suite = unittest.makeSuite(test_SBMLFunctions)
else:
    message = no_libsbml_msg

if __name__ == '__main__':
    if _HAVE_SBML:
        unittest.main()
    else:
        print no_libsbml_msg
