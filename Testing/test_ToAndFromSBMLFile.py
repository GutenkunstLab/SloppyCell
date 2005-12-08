import copy
import os
import unittest

from SloppyCell.ReactionNetworks import *
# Check whether we actually have the SBML methods.
_HAVE_SBML = (hasattr(IO, 'to_SBML_file') and hasattr(IO, 'from_SBML_file'))

from TestNetwork import net
net = copy.deepcopy(net)
net.compile()
 
class test_SBMLFunctions(unittest.TestCase):
    def test_ToAndFromSBMLFile(self):
        """Testing the reading and writing of an SBML file from and to a network"""
        filename = 'net.sbml'
	IO.to_SBML_file(net, filename)
	netread = IO.from_SBML_file(filename)
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
