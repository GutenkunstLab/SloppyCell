import os
import scipy
import unittest

from SloppyCell.ReactionNetworks import *

base_net = Network_mod.Network('test')
base_net.add_compartment('testTube')

base_net.add_species('x', 'testTube', 1000.0)
base_net.add_species('y', 'testTube', 0.0)
base_net.add_parameter('A')
base_net.addReaction('x->y', kineticLaw='A*x', stoichiometry={'x':-1, 'y':+1})

base_net.add_parameter('B')
base_net.addReaction('y->x', kineticLaw='B*y', stoichiometry={'y':-1, 'x':+1})

times = scipy.arange(0., 5., .1)
params = [1., 1.]

class test_Events(unittest.TestCase):
    def test_basic(self):
        """Test that the stochastic integrator's basic functions work"""
        net = base_net.copy('test_basic')
        net.set_deterministic()
        resDt = net.calculate({'x':times}, params)

        net.set_stochastic()
        resSt = net.calculate({'x':times}, params)

        import pylab
        pylab.figure()
        pylab.plot(times, [resDt['x'][_] for _ in times], 'b--')
        pylab.plot(times, [resSt['x'][_] for _ in times], 'b-')
        
        pylab.xlabel('Time')
        pylab.ylabel('Number of Molecules')
        pylab.title('Stochastic Test - Do the trajectories look similar?')

    def test_RngSeed(self):
        """Test that the reseeding the same RNG seed produces identical results"""
        net = base_net.copy('test_RngSeed')
        net.set_stochastic(seed=24048296)

        res1 = net.calculate({'x':times}, params)

        net.stochastic['reseed']=True
        res2 = net.calculate({'x':times}, params)

        for t in times:
            self.assertEqual(res1['x'][t], res2['x'][t],
                             'Same RNG seed failed to produce the same result')

    def test_Assignment(self):
        """Test that assignments work in the stochastic integrator"""
        net = base_net.copy('test_Assignment')

        net.add_species('u', 'testTube')
        net.add_assignment_rule('u', 'y*x+sqrt(time)')

        net.add_species('v', 'testTube', 1000.0)
        net.addReaction('v degd', kineticLaw='u*v', stoichiometry={'v':-1})
        
        net.set_stochastic()

        res = net.calculate({'x':times, 'y':times, 'u':times, 'v':times},
                            params)

        total = res['x'][0.0]+res['y'][0.0]
        resTimes = res['x'].keys()
        resTimes.sort()
        for t in resTimes:
            self.assertAlmostEqual(res['x'][t]+res['y'][t], total, 6,
                                   'Failed to conserve mass.')
            for var in ['x','y','u','v']:
                self.assertTrue(res[var][t]>=0,
                                'Failed to keep %s concentration above 0.'%var)

    def test_NoReactions(self):
        """Test that the stochastic integrator is quiet when no reactions are found"""
        net = Network_mod.Network('test')
        net.add_compartment('testTube')

        net.add_species('x', 'testTube', 1000.0)
        net.add_species('y', 'testTube', 0.0)

        net.add_parameter('A')
        net.add_parameter('B')

        # These wouldn't be used (no rate rules are), but do raise warnings
        #net.add_rate_rule('x', '-A*x+B*y')
        #net.add_rate_rule('y', 'A*x-B*y')

        net.set_stochastic()

        res = net.calculate({'x': times}, params)

        for time, val in res['x'].items():
            self.assertTrue(val==1000.0,
                            'Failed to silently fail with no rxns present.')
        
suite = unittest.makeSuite(test_Events)
if __name__ == '__main__':
    unittest.main()
