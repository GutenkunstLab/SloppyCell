import unittest

from SloppyCell.ReactionNetworks import *
from TestNetwork import net, m, params

m.cost(params)

class test_Plotting(unittest.TestCase):
    def test_data_plot(self):
        """Ensure plotting basic data doesn't crash"""
        Plotting.plot_model_data(m)

    def test_results_plot(self):
        """Ensure plotting basic results doesn't crash"""
        Plotting.plot_model_results(m)

suite = unittest.makeSuite(test_Plotting, 'test')

if __name__ == '__main__':
    unittest.main()
