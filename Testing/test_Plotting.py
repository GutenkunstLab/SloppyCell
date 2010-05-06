import unittest

from SloppyCell.ReactionNetworks import *
from TestNetwork import net, m, params

class test_Plotting(unittest.TestCase):
    def test_data_plot(self):
        """Ensure plotting basic data doesn't crash"""
        m.cost(params)
        Plotting.plot_model_data(m)

    def test_results_plot(self):
        """Ensure plotting basic results doesn't crash"""
        m.cost(params)
        Plotting.plot_model_results(m)

# Only run if Plotting available
import SloppyCell
_HAVE_PLOTTING = hasattr(SloppyCell, 'Plotting')
no_plotting_msg = 'Plotting not available.'
if _HAVE_PLOTTING:
    suite = unittest.makeSuite(test_Plotting)
else:
    message = no_plotting_msg
    
if __name__ == '__main__':
    if _HAVE_PLOTTING:
        unittest.main()
    else: 
        print no_plotting_msg

