import unittest

import scipy

from SloppyCell.ReactionNetworks import *

base_net = Network('test')
base_net.add_compartment('blank')
base_net.add_species('X', 'blank', 1)
base_net.add_species('Y', 'blank', 0.5)

class test_FreeEnergy(unittest.TestCase):
    def test_uniform_free_sf(self):
        """Test that single floating scale factor handled correctly"""
        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (5.0, 0.5),
                                      2.0: (5.0, 1.0)}}})
        m = Model([expt], [net])

        c = m.cost([])
        self.assertEqual(c, 0, 'Failed on zero-cost integration.')

        # ak is sum(y**2/sigma**2)
        ak = 1.0**2/0.5**2 + 1.0**2/1.0**2

        for T in [1, 0.137]:
            entropy = scipy.log(scipy.sqrt(2*scipy.pi*T/ak))
            self.assertAlmostEqual(m.free_energy([], T), c-T*entropy, 4,
                                   'Failed on T=%f free energy' % T)

    def test_fixed_sf(self):
        """Test that single fixed scale factor handled correctly"""
        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (4.0, 0.5),
                                      2.0: (5.0, 1.0)}}})
        expt.set_fixed_sf({'X': 3.0})
        m = Model([expt], [net])

        c = 0.5*( ((4.0 - 3.0)/0.5)**2 + ((5.0 - 3.0)/1.0)**2)
        self.assertAlmostEqual(m.cost([]), c, 3, 'Failed on cost integration.')

        for T in [1, 0.137]:
            entropy = 0
            self.assertAlmostEqual(m.free_energy([], T), c-T*entropy, 4,
                                   'Failed on T=%f free energy' % T)

    def test_shared_sf(self):
        """Test that shared floating scale factor handled correctly"""
        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (4.0, 0.5),
                                      2.0: (5.0, 1.0)},
                                'Y': {0.5: (3.0, 0.25)}}})
        expt.set_shared_sf([('X', 'Y')])
        m = Model([expt], [net])

        # ak is sum(y**2/sigma**2)
        ak = 1.0**2/0.5**2 + 1.0**2/1.0**2 + 0.5**2/0.25**2
        # bk is sum(y*d/sigma**2)
        bk = 1.0*4.0/0.5**2 + 1.0*5.0/1.0**2 + 0.5*3.0/0.25**2
        B = bk/ak

        c = 0.5*(((B*1.0 - 4.0)/0.5)**2 + ((B*1.0 - 5.0)/1.0)**2
                 + ((B*0.5 - 3.0)/0.25)**2)
        self.assertAlmostEqual(m.cost([]), 4.0, 4, 'Failed on shared sf cost '
                               'calc.')

        for T in [1, 0.137]:
            entropy = scipy.log(scipy.sqrt(2*scipy.pi*T/ak))
            self.assertAlmostEqual(m.free_energy([], T), c-T*entropy, 4,
                                   'Failed on T=%f free energy' % T)

    def test_log_gaussian_prior_sf(self):
        """
        Test that single scale factor with a log gaussian prior is handled
        correctly
        """
        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (4.0, 0.5),
                                      2.0: (5.0, 1.0)}}})
        expt.set_sf_prior('X', 'gaussian in log sf', 
                          (scipy.log(2.0), scipy.log(1e3)))
        m = Model([expt], [net])

        # ak is sum(y**2/sigma**2)
        ak = 1.0**2/0.5**2 + 1.0**2/1.0**2
        # bk is sum(y*d/sigma**2)
        bk = 1.0*4.0/0.5**2 + 1.0*5.0/1.0**2

        c = m.cost([])

        # The comparison entropies were calculated in Mathematica.
        T = 1.0
        entropy = (c - m.free_energy([],T))/T
        self.assertAlmostEqual(entropy, -1.31481, 4,
                               'Failed on T = 1 entropy')

        T = 0.137
        entropy = (c - m.free_energy([],T))/T
        self.assertAlmostEqual(entropy, -2.31894, 4,
                               'Failed on T = 0.137 entropy')

    def test_bad_prior_type(self):
        """
        Test that an unknown prior type raises and exception.
        """
        T = 0.876

        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (4.0, 0.5),
                                      2.0: (5.0, 1.0)},
                                'Y': {0.5: (3.0, 0.25)}}})
        m = Model([expt], [net])

        self.assertRaises(ValueError, expt.set_sf_prior, 
                          'X', 'bad prior type')

    def test_changing_sf_priors(self):
        """
        Test combinations of log and uniform sf priors.
        """
        T = 0.876

        net = base_net.copy()
        expt = Experiment()
        expt.set_data({'test': {'X': {1.0: (4.0, 0.5),
                                      2.0: (5.0, 1.0)},
                                'Y': {0.5: (3.0, 0.25)}}})
        m = Model([expt], [net])

        # First we calculate the entropies for various individual scale
        # factors. We do this by fixing the other scale factor.
        expt.set_fixed_sf({'Y': 1.0})
        ent_X_uniform = (m.cost([]) - m.free_energy([], T))/T

        expt.set_fixed_sf({'X': 1.0})
        ent_Y_uniform = (m.cost([]) - m.free_energy([], T))/T

        expt.set_sf_prior('X', 'gaussian in log sf', 
                          (scipy.log(2.0), scipy.log(1e3)))
        expt.set_fixed_sf({'Y': 1.0})
        ent_X_log = (m.cost([]) - m.free_energy([], T))/T

        expt.set_sf_prior('Y', 'gaussian in log sf', 
                          (scipy.log(2.0), scipy.log(1e3)))
        expt.set_fixed_sf({'X': 1.0})
        ent_Y_log = (m.cost([]) - m.free_energy([], T))/T

        # Now we try various combinations.
        expt.set_fixed_sf({})
        ent_both_log = (m.cost([]) - m.free_energy([], T))/T
        self.assertAlmostEqual(ent_both_log, ent_X_log+ent_Y_log, 6)

        expt.set_sf_prior('Y', 'uniform in sf', None)
        ent_this = (m.cost([]) - m.free_energy([], T))/T
        self.assertAlmostEqual(ent_this, ent_X_log+ent_Y_uniform, 6)

        expt.set_sf_prior('X', 'uniform in sf', None)
        ent_this = (m.cost([]) - m.free_energy([], T))/T
        self.assertAlmostEqual(ent_this, ent_X_uniform+ent_Y_uniform, 6)

suite = unittest.makeSuite(test_FreeEnergy)

if __name__ == '__main__':
    unittest.main()
