import copy
import unittest

import scipy

from SloppyCell.ReactionNetworks import *
import SloppyCell
import TestNetwork
import SloppyCell.ReactionNetworks.OptDesign as OD
import os
import SloppyCell.Utility
save = SloppyCell.Utility.save
load = SloppyCell.Utility.load

net = copy.deepcopy(TestNetwork.net)
m = copy.deepcopy(TestNetwork.m)
p = m.params.__copy__()

class test_OptDesign(unittest.TestCase):
    def test_OptDesign_funcs(self):
        """Test that the functions in OptDesign work correctly"""
        save(p,'TNparams')
        j,jtj = m.GetJandJtJInLogParameters(p)
        save(jtj,'TNjtj')
        # Make the nensitivity trajectory for this network and save
        # it in 'TNsenstraj'
        times = scipy.linspace(0, 2.6, 100)
        times = scipy.append(times, 1.1)
        times.sort()
        OD.make_sens_traj(net,p,times,'TNsenstraj')
        OD.setup('TNparams',net,'TNsenstraj','TNjtj')
        best_change,best_chem,best_time = OD.design_over_chems(['x'],['y'],
                scipy.log(1000.0))
        times, bestfit, var = OD.variances(['x','y'],scipy.log(1000.0))
        times, bestfit, var = OD.variances_log_chems(['x','y'],scipy.log(1000.0))
        sens_vect = OD.get_sens_vect('x',1.1)
        sens_array = OD.get_sens_array('x')
        # tidy up
        os.unlink('TNjtj')
        os.unlink('TNparams')
        os.unlink('TNsenstraj')

suite = unittest.makeSuite(test_OptDesign)

if __name__ == '__main__':
    unittest.main()
