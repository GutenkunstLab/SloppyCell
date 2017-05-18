import unittest
import copy
import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net
tlist = scipy.array([0] + [0.8*x for x in range(1, 51)])

class test_IntegrateWithoutDerivatives(unittest.TestCase):
    def test_basic(self):
        local_net = copy.deepcopy(algebraic_net)

        # need to add a new parameter to change the network structure and
        # force it to recompile
        local_net.addParameter(id='dummy_par',value=1.0)
        local_net.disable_deriv_funcs()
        local_net.compile()

        funcs_no_derivs = ['res_function', 'alg_deriv_func', 'alg_res_func',\
                           'integrate_stochastic_tidbit', 'root_func']
        self.assertEqual(local_net._dynamic_funcs_python.keys(),
                         funcs_no_derivs)
        traj = Dynamics.integrate(local_net, tlist)

        self.assertAlmostEqual(traj.get_var_val('X0',4.8), 
                               0.618783392, 5)
        self.assertAlmostEqual(traj.get_var_val('X1',21.6), 
                               0.653837775, 5)
        
suite = unittest.makeSuite(test_IntegrateWithoutDerivatives)
if __name__ == '__main__':
    unittest.main()
