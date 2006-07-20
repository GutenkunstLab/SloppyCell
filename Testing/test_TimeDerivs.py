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


class test_TimeDerivs(unittest.TestCase):
    def test_TimeDerivs(self):
        """Test time derivatives of the trajectory and the sensitivity trajectory """
        net.add_int_times = False
        # need really finely spaced time points for the finite differencing to work
        Traj = Dynamics.integrate(net,scipy.linspace(0.0,.001,1000), return_derivs=True)
        sensTraj = Dynamics.integrate_sensitivity(net,scipy.linspace(0.0,.001,1000), return_derivs=True)
        x_vals = Traj.getVariableTrajectory('x')
        x_vals_deriv = scipy.diff(x_vals)/scipy.diff(Traj.timepoints)
        relerror = scipy.linalg.norm( x_vals_deriv-Traj.getVariableTrajectory(('x','time'))[:-1] )/scipy.linalg.norm(x_vals_deriv)

        x_vals_wrtA = sensTraj.getVariableTrajectory(('x','A'))
        x_vals_wrtA_deriv = scipy.diff(x_vals_wrtA)/scipy.diff(sensTraj.timepoints)

        relerror2 = scipy.linalg.norm( x_vals_wrtA_deriv-sensTraj.getVariableTrajectory(('x','A','time'))[:-1] )\
                    /scipy.linalg.norm(x_vals_wrtA_deriv)

        self.assertAlmostEqual(relerror, 0.0, 5,'Failed on time deriv of x in Traj')
        self.assertAlmostEqual(relerror2, 0.0, 5,'Failed on time deriv of sensitivity of x wrt A in sensTraj')

suite = unittest.makeSuite(test_TimeDerivs)

if __name__ == '__main__':
    unittest.main()
