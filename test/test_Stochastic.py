import os, sys
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

base_net.add_event('event 0', 'gt(time, 5)', {'x': 0.0, 'y': 1000.0})

times = scipy.arange(0., 10.05, .1)
params = [1., 1.]

class test_Events(unittest.TestCase):
    def test_basic(self):
        """Test that the stochastic integrator's basic functions work"""
        net = base_net.copy('test_basic')
        net.set_deterministic()
        net.Calculate({'x':times}, params)
        trajDt = net.trajectory.time_slice(0.,10.)
        
        net.set_stochastic()
        net.Calculate({'x':times}, params)
        trajSt = net.trajectory.time_slice(0.,10.)
        
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

    def test_NonBuildingStochastic(self):
        """Test that the stochastic integrator does not build when unsupported network components are present"""
        def passed(_net, msg):
            _net.compile()
            self.assertTrue(
                'pass' in 
                _net._dynamic_funcs_python.get('integrate_stochastic_tidbit'),
                'Stochastic integrator compiled with %s.'%msg)
        
        net = base_net.copy('test_Network')

        # No Reactions
        for key in net.reactions.keys():
            net.remove_component(key)
        passed(net, 'no reactions')
        net.addReaction('y->0', kineticLaw='-B*y', stoichiometry={'y':-1})

        # Rate Rules
        net.add_rate_rule('x', '-A*x')
        passed(net, 'a rate rule')
        net.rateRules.remove_by_key('x')

        # Events
        net.add_event('event 1', 'gt(y, 5)', {'A': 2})
        passed(net, 'an event')
        net.remove_component('event 1')
        
        # Algebraic Rules / Vars
        net.add_algebraic_rule('x+y')
        passed(net, 'an algebraic rule')
        net.algebraicRules.remove_by_key('x+y')

        # Non-int-able stoichiometry
        net.addReaction('x->0', kineticLaw='-A*x', stoichiometry={'x':1.001})
        passed(net, 'non-integer stoichiometry')
        net.remove_component('x->0')

        net.addReaction('x->0', kineticLaw='-A*x', stoichiometry={'x':'2*x'})
        passed(net, 'a mathematical stoichiometric expression')

        # Finally, test that the network fails predictably when integrated
        net.set_stochastic()
        err = False
        try:
            net.calculate({'x':times}, params)
        except RuntimeError:
            self.assertTrue("# Noncastable stoichiometry for x->0: '2*x'" in
                            '%s'%sys.exc_value,
                            'Raised an unknown RuntimeError: %s'%sys.exc_value)
            err = True
        self.assertTrue(err, 'Failed to raise error with unsupported reaction')
        
suite = unittest.makeSuite(test_Events)
if __name__ == '__main__':
    unittest.main()
