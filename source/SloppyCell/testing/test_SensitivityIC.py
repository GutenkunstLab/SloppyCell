import unittest

import scipy,copy
from SloppyCell.ReactionNetworks import *

from TestNetwork import net
net = copy.deepcopy(net)
net.setInitialVariableValue('A', 1.0)
net.setInitialVariableValue('B', 2.0)
net.addParameter('xIC',4.0)
net.setInitialVariableValue('x','2.0*exp(xIC)')
net.compile()


class test_SensitivityIC(unittest.TestCase):
    def test_SensitivityIC(self):
        """Test that sensitivity vector is initialized correctly"""
	sensTraj = net.integrateSensitivity(scipy.linspace(0.0,1.0,10))
	x_wrt_A = sensTraj.getVariableTrajectory(('x','A'))[0]
	x_wrt_B = sensTraj.getVariableTrajectory(('x','B'))[0]
	x_wrt_xIC = sensTraj.getVariableTrajectory(('x','xIC'))[0]
	self.assertAlmostEqual(x_wrt_A, 0.0, 8,'Failed on x_wrt_A sensitivity at t=0')
	self.assertAlmostEqual(x_wrt_B, 0.0, 8,'Failed on x_wrt_B sensitivity at t=0')
	self.assertAlmostEqual(x_wrt_xIC, 2.0*scipy.exp(4.0), 8,'Failed on x_wrt_xIC sensitivity at t=0')


suite = unittest.makeSuite(test_SensitivityIC)

if __name__ == '__main__':
    unittest.main()
