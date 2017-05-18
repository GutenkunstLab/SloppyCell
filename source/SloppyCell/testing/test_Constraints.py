import unittest
import os

import scipy

from SloppyCell.ReactionNetworks import *

from AlgTestNets import algebraic_net_constraints
base_net = algebraic_net_constraints.copy()
tlist_algebraic_net = scipy.array([0] + [0.8*x for x in range(1, 51)])

class test_Constraints(unittest.TestCase):
    def test_contraint_fire(self):
        """
        Test that a constraint violation causes an exception.
        """
        
        net = base_net.copy('test')

        self.assertRaises((Utility.ConstraintViolatedException),
                          Dynamics.integrate, net, tlist_algebraic_net,
                          fill_traj=True, redirect_msgs=True)


    def test_constraint_time(self):
        """
        Test that the constraint violation happens at the right time.
        """        
        net = base_net.copy('test')
        try:
            traj = Dynamics.integrate(net, tlist_algebraic_net)
        except Utility.ConstraintViolatedException, cve:
            self.assertAlmostEqual(cve.time, 16.436798814)
            self.assertEqual(cve.message, 'X1 is big!')            

    def test_add_constraint(self):
        """
        Test that multiple constraints work.
        """          
        net = base_net.copy('test')

        net.addConstraint('TimeTooLarge','lt(time,1.0)','Time got too large')

        self.assertRaises((Utility.ConstraintViolatedException),
                          Dynamics.integrate, net, tlist_algebraic_net,
                          fill_traj=True, redirect_msgs=True)

        try:
            traj = Dynamics.integrate(net, tlist_algebraic_net)
        except Utility.ConstraintViolatedException, cve:
            self.assertAlmostEqual(cve.time, 1.0)

    def test_constraints_off(self):
        """
        Test that events fire properly if the constraints are turned off.
        """
        net = base_net.copy('test')

        net.addConstraint('TimeTooLarge','lt(time,1.0)','Time got too large')

        traj = Dynamics.integrate(net, tlist_algebraic_net, use_constraints=False)

        event_indeces = [net.events.index_by_key('event1'),net.events.index_by_key('event2')]

        (te,ye,ie) = traj.event_info

        for e_ind in event_indeces:
            self.assertEqual(True, e_ind in ie)
        for t1, t2 in zip(te, [13.30131485,
                              18.0]):
            self.assertAlmostEqual(t1, t2)

    def test_constraint_chaining(self):
        """
        Test that an event firing can cause a constraint to fire (by assigning
        model variables such that the model is in an invalid state).
        """
        net = base_net.copy('test')

        eventAssignments = {'X1':10}

        net.addEvent('increase_X1', 'geq(time,2.0)', eventAssignments)

        self.assertRaises((Utility.ConstraintViolatedException),
                          Dynamics.integrate, net, tlist_algebraic_net,
                          fill_traj=True, redirect_msgs=True)

        try:
            traj = Dynamics.integrate(net, tlist_algebraic_net)
        except Utility.ConstraintViolatedException, cve:
            self.assertAlmostEqual(cve.time, 2.0)

    def test_constraints_time_zero(self):
        """
        If checking for constraint violations at time zero.
        """
        net = base_net.copy('test')

        # set X1 to be too big for the constraint
        net.set_var_ic('X1', 0.6)

        self.assertRaises((Utility.ConstraintViolatedException),
                          Dynamics.integrate, net, tlist_algebraic_net,
                          fill_traj=True, redirect_msgs=True)

        try:
            traj = Dynamics.integrate(net, tlist_algebraic_net)
        except Utility.ConstraintViolatedException, cve:
            self.assertAlmostEqual(cve.time, 0.0)

        
suite = unittest.makeSuite(test_Constraints)
if __name__ == '__main__':
    unittest.main()
