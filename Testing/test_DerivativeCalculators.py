import copy
import unittest

import scipy

from SloppyCell.ReactionNetworks import *
import SloppyCell
import TestNetwork

net = TestNetwork.net
m = copy.deepcopy(TestNetwork.m)
# the following is useful to see if prior derivatives are computed
# correctly
m.AddResidual(Residuals.PriorInLog('priorOnA', 'A',
                                   scipy.log(1.1),scipy.log(2.0) ) )

class test_DerivativeCalculators(unittest.TestCase):
    def test_jacobian_sens(self):
        """
        Test that sensitivity and finite-difference jacobians agree.
        """
        params = m.get_params()

        J_sens = m.jacobian_sens(params)
        J_fd = m.jacobian_fd(params, 1e-4)

        for res_name in J_sens.keys():
            sens_vals = J_sens.get(res_name)
            fd_vals = J_sens.get(res_name)
            for sens_val, fd_val in zip(sens_vals, fd_vals):
                if (sens_val == 0) and (fd_val == 0):
                    continue
                rel_diff = abs(sens_val - fd_val)/max(abs(sens_val),abs(fd_val))
                self.assertAlmostEqual(rel_diff, 0, 3, 
                                       'Failed on residual %s.' % str(res_name))

    def test_gradient(self):
        """
        Test sensitivity gradient of cost.
        """
        params = m.get_params()

        grad_sens = m.gradient_sens(params)

        c0 = m.cost(params)
        eps = 1e-4
        for p_key, p_val in params.items():
            pplus = params.copy()
            pplus.set(p_key, p_val + eps)
            cplus = m.cost(pplus)
            deriv = (cplus - c0)/eps
            self.assertAlmostEqual(deriv, grad_sens.get(p_key), 2)

    def test_gradient_log_params(self):
        """
        Test sensitivity gradient wrt log params of cost.
        """
        params = m.get_params()

        grad_sens = m.gradient_log_params_sens(scipy.log(params))

        c0 = m.cost(params)
        eps = 1e-6
        for p_key, p_val in params.items():
            pplus = params.copy()
            pplus.set(p_key, scipy.exp(scipy.log(p_val) + eps))
            cplus = m.cost(pplus)
            deriv = (cplus - c0)/eps
            self.assertAlmostEqual(deriv, grad_sens.get(p_key), 2)
	 
suite = unittest.makeSuite(test_DerivativeCalculators)

if __name__ == '__main__':
    unittest.main()
