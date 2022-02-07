import unittest 

import scipy
import SloppyCell.Utility as Utility
import SloppyCell.daskr
from SloppyCell.daskr import daeint

redir = Utility.Redirector()

################################################################################

# Van Der Pol oscillator equations
# This test problem is from the daskr documentation.
def vdp_func(y, t):
    ydot = scipy.zeros(2, scipy.float_)
    ydot[0] = y[1]
    ydot[1] = 100*(1-y[0]**2)*y[1] - y[0]
    return ydot
    
def vdp_res_func(t, y, yprime, rpar):
    return vdp_func(y, t) - yprime 


def vdp_Dfun(t, y, yprime, cj, rpar):
    pd = scipy.zeros((2,2), scipy.float_)
    pd[0,0] = -cj
    pd[0,1] = 1
    pd[1,0] = -2*100*y[0]*y[1]-1
    pd[1,1] = 100*(1-y[0]**2)-cj
    return pd

def vdp_rt_func(t, y, yp, rpar):
    trigger = y[0]
    return scipy.asarray([trigger])

# initial settings for vdp system
y0_vdp = scipy.array([2.0, 0])
tlist_vdp = scipy.array([0] + [20*x for x in range(1, 11)])
t0_vdp = tlist_vdp[0]
yp0_vdp = vdp_func(y0_vdp, t0_vdp)
num_events_vdp = 1
abstol_vdp = scipy.array([0.1e-5, 0.1e-3])
reltol_vdp = scipy.array([0.1e-5, 0.1e-5])

################################################################################

# AlgebraicRules_BasicExample from the SBML Level 2 Version 1 test suite
# Variables x, y

def algExampleBasic_func(y,t):
    ydot = scipy.zeros(3, scipy.float_)
    ydot[0] = 1*y[1]
    return ydot

def algExampleBasic_res_func(t, y, yprime, rpar):
    res = scipy.zeros(3, scipy.float_)
    ypcalc = algExampleBasic_func(y,t)
    res[0] = ypcalc[0] - yprime[0]
    res[1] = -y[2]+y[0]+y[1]
    res[2] = yprime[2]         
    return res

# initial settings for basic algebraic system
y0_algBasic = scipy.array([0.5, 0.5, 1])
tlist_algBasic = scipy.array([0] + [0.2*x for x in range(1, 51)])
t0_algBasic = tlist_algBasic[0]
yp0_algBasic = algExampleBasic_func(y0_algBasic, t0_algBasic)

abstol_algBasic = scipy.array([0.1e-8, 0.1e-8, 0.1e-8])
reltol_algBasic = scipy.array([0.1e-5, 0.1e-5, 0.1e-5])


################################################################################

# AlgebraicRules_FastReactionExample from the SBML Level 2 Version 1 test suite
# The given assignmentRule is made into an algebraic rule
# Variables X0, X1, T, S1, S2

#Parameters
Keq = 2.5
k1 = 0.1
k2 = 0.15

def algExample_func(y,t):
    ydot = scipy.zeros(5, scipy.float_)
    ydot[0] = -k1*y[0]
    ydot[1] = k2*y[4]
    ydot[2] = k1*y[0] - k2*y[4]
    return ydot

def algExample_res_func(t, y, yprime, ires):
    res = scipy.zeros(5, scipy.float_)
    ypcalc = algExample_func(y,t)
    res[0] = ypcalc[0] - yprime[0]
    res[1] = ypcalc[1] - yprime[1]
    res[2] = ypcalc[2] - yprime[2]
    res[3] = (y[3] + y[4] - y[2])
    res[4] = (y[4] - Keq*y[3])
    return res

# initial settings for algebraic fast reaction system
y0_alg = scipy.array([1.0, 0, 0, 0, 0])
tlist_alg = scipy.array([0] + [0.8*x for x in range(1, 51)])
t0_alg = tlist_alg[0]
yp0_alg = algExample_func(y0_alg, t0_alg)

num_events = 1
abstol_alg = scipy.array([0.1e-5, 0.1e-5, 0.1e-5, 0.1e-5, 0.1e-5])
reltol_alg = scipy.array([0.1e-5, 0.1e-5, 0.1e-5, 0.1e-5, 0.1e-5])


################################################################################

# Simple linear equation
# This test problem is for testing tstop
# Note:  Some time points in the tlist will be skipped when tstop is
# encountered.

def linear_func(y, t):
    ydot = scipy.zeros(1, scipy.float_)
    ydot[0] = -5
    return ydot
    
def linear_res_func(t, y, yprime, ires):
    return linear_func(y, t) - yprime 

# initial settings for simple linear system
y0_linear = scipy.array([100])
tlist_linear = scipy.array([0] + [20*x for x in range(1, 11)])
t0_linear = tlist_linear[0]
yp0_linear = linear_func(y0_linear, t0_linear)
abstol_linear = scipy.array([0.1e-5])
reltol_linear = scipy.array([0.1e-5])
tstop_linear = 201


################################################################################

# The non_neg example tests the checking of y going negative
# This system has a rapid change in dynamics at y = k, and it's
# easy for the integratory to miss if non-negativity is not enforced.

#Parameters
k = 1e-12

def non_neg_func(y,t):
    ydot = scipy.zeros(1, scipy.float_)
    ydot[0] = -y[0]/(k+y[0])
    return ydot

def non_neg_res_func(t, y, yprime, ires):
    res = scipy.zeros(1, scipy.float_)
    ypcalc = non_neg_func(y,t)
    res[0] = ypcalc[0] - yprime[0]
    return res


# initial settings for basic non negative system
y0_non_neg = scipy.array([1.0])
tlist_non_neg = scipy.array([0] + [0.04*x for x in range(1, 51)])
t0_non_neg = tlist_non_neg[0]
yp0_non_neg = non_neg_func(y0_non_neg, t0_non_neg)

abstol_non_neg = scipy.array([0.1e-5])
reltol_non_neg = scipy.array([0.1e-5])

################################################################################

# Simple time dependent trigonometric system
# This test problem is for testing tstop
# Note:  Some time points in the tlist will be skipped when tstop is
# encountered.

def trig_func(y, t):
    ydot = scipy.zeros(1, scipy.float_)
    ydot[0] = scipy.cos(t)
    return ydot
    
def trig_res_func(t, y, yprime, ires):
    return trig_func(y, t) - yprime 

# initial settings for simple linear system
y0_trig = scipy.array([0])
tlist_trig = scipy.array([0] + [1000*x for x in range(1, 3)])
t0_trig = tlist_trig[0]
yp0_trig = trig_func(y0_linear, t0_linear)
abstol_trig = scipy.array([0.1e-5])
reltol_trig = scipy.array([0.1e-5])


################################################################################


class test_daskr(unittest.TestCase):
    def test_basic(self):
        """ Basic test of daskr """
        y, t, ypout, t_root, y_root, i_root = daeint(vdp_res_func, tlist_vdp,
                                                     y0_vdp, yp0_vdp,
                                                     rtol = reltol_vdp,
                                                     atol = abstol_vdp,
                                                     intermediate_output=False)

        self.assertAlmostEqual(y[1][0], 1.85821444, 4)
        self.assertAlmostEqual(y[3][0], 0.1484599E+01, 4)
        self.assertAlmostEqual(y[7][0], -0.1501730E+01, 4)
        self.assertAlmostEqual(y[10][0], 0.1718428E+01, 4)

        self.assertAlmostEqual(y[2][1], -0.9068522E-02, 3)
        self.assertAlmostEqual(y[4][1], -0.5847012E-01, 3)
        self.assertAlmostEqual(y[8][1], 0.3569131E-01, 3)
        self.assertAlmostEqual(y[9][1], -0.7422161E-02, 3)

    def test_Dfun(self):
        """ Test user-supplied Jacobian """
        y, t, ypout, t_root, y_root, i_root = daeint(vdp_res_func, tlist_vdp,
                                                     y0_vdp, yp0_vdp,
                                                     jac=vdp_Dfun,
                                                     rtol = reltol_vdp,
                                                     atol = abstol_vdp,
                                                     intermediate_output=False)
        
        self.assertAlmostEqual(y[1][0], 1.85821444, 4)
        self.assertAlmostEqual(y[6][1], 8.93022e-3, 4)

    def test_term_roots(self):
        """ Test root finding with termination """
        y, t, ypout, t_root, y_root, i_root = daeint(vdp_res_func, tlist_vdp,
                                                     y0_vdp, yp0_vdp,
                                                     nrt=1, 
                                                     rt=vdp_rt_func,
                                                     rtol = reltol_vdp,
                                                     atol = abstol_vdp,
                                                     intermediate_output=False)
  
        self.assertAlmostEqual(t_root, 0.8116351E+02, 4)
        self.assertAlmostEqual(y_root[0], -0.3295063E-12, 4)
        self.assertAlmostEqual(y_root[1], -0.6714100E+02, 3)
        self.assertEqual(i_root[0], -1)

    def test_tstop(self):
        """ Test that integration will not continue past tstop """
        y, t, ypout, t_root, y_root, i_root = daeint(linear_res_func, 
                                                     tlist_linear,
                                                     y0_linear, yp0_linear,
                                                     rtol=reltol_linear,
                                                     atol=abstol_linear,
                                                     tstop=tstop_linear)

        
        # Check that the final time point returned is for tstop
        self.assertAlmostEqual(t[-1], tstop_linear, 4)
        self.assertAlmostEqual(y[2][0], -100, 4)


    
    def test_algebraic_basic(self):
        """ Test a simpler dae system (algebraicRules-basic-l2.xml) """
        y, t, ypout, t_root, y_root, i_root = daeint(algExampleBasic_res_func, 
                                                     tlist_algBasic,
                                                     y0_algBasic, yp0_algBasic,
                                                     rtol = reltol_algBasic,
                                                     atol = abstol_algBasic)

        self.assertAlmostEqual(y[1][0], 0.590635382065755, 4)
        self.assertAlmostEqual(y[13][0], 0.962863096631099, 4)

        self.assertAlmostEqual(y[15][1], 0.0248936510867585, 4)
        self.assertAlmostEqual(y[27][1], 0.00225832507503575, 4)



    def test_algebraic_fastreactionexample(self):
        """ Test a dae system (algebraicRules-fastReactionExample-l2.xml) """
        y, t, ypout, t_root, y_root, i_root = daeint(algExample_res_func, 
                                                     tlist_alg,
                                                     y0_alg, yp0_alg,
                                                     rtol = reltol_alg,
                                                     atol = abstol_alg)

        self.assertAlmostEqual(y[1][0], 0.9231163463, 4)
        self.assertAlmostEqual(y[13][0], 0.353454681, 4)

        self.assertAlmostEqual(y[8][1], 0.142837751, 4)
        self.assertAlmostEqual(y[20][1], 0.492844600, 4)

        self.assertAlmostEqual(y[15][2], 0.346376313, 4)
        self.assertAlmostEqual(y[27][2], 0.230837103, 4)
        
        self.assertAlmostEqual(y[22][3], 0.081296859, 4)
        self.assertAlmostEqual(y[37][3], 0.039501126, 4)

        self.assertAlmostEqual(y[29][4], 0.150075280, 4)
        self.assertAlmostEqual(y[41][4], 0.078591978, 4)

        self.assertAlmostEqual(y[50][0], 0.018315639, 4)
        self.assertAlmostEqual(y[50][1], 0.917958431, 4)
        self.assertAlmostEqual(y[50][2], 0.06372593, 4)
        self.assertAlmostEqual(y[50][3], 0.018207409, 4)
        self.assertAlmostEqual(y[50][4], 0.045518522, 4)


    def test_maxsteps_on(self):
        """ Test to make sure the max_steps parameter works """
        y, t, ypout, t_root, y_root, i_root = daeint(trig_res_func, tlist_trig,
                                                     y0_trig, yp0_trig,
                                                     rtol = reltol_trig,
                                                     atol = abstol_trig,
                                                     max_steps = 7500)

        # the integrator will only get to the specified time points if
        # max_steps is increased significantly above the default
        self.assertAlmostEqual(y[1][0], 0.82689894, 4)
        self.assertAlmostEqual(y[2][0], 0.93004774, 4)

    def test_maxsteps_off(self):
        """ Test to make sure the trig_func problem will cause an error \
if max_steps is not set """

        redir.start()
        try:
            self.assertRaises(SloppyCell.daskr.daeintException, 
                              daeint(trig_res_func, tlist_trig,
                                     y0_trig, yp0_trig,
                                     rtol = reltol_trig,
                                     atol = abstol_trig))
        except SloppyCell.daskr.daeintException:
            pass
        messages = redir.stop()



    def test_algebraic_calculate_ic(self):
        """ Test automatic calculation of initial conditions """

        # pass an inconsistent set of initial conditions to the fast reaction
        # example
        y0_inconsistent = scipy.array([1.0, 0, 0, 1500, 15])
        yp0_inconsistent = algExample_func(y0_inconsistent, t0_alg)
        var_types_inconsistent = scipy.array([1, 1, 1, -1, -1])
        
        y, t, ypout, t_root, y_root, i_root = daeint(algExample_res_func, 
                                                     tlist_alg,
                                                     y0_inconsistent, yp0_alg,
                                                     rtol = reltol_alg,
                                                     atol = abstol_alg,
                                                     calculate_ic = True,
                                                     var_types = var_types_inconsistent)

        # check to make sure the initial condition was calculated correctly
        self.assertAlmostEqual(y[0][0], 1., 4)
        self.assertAlmostEqual(y[0][1], 0., 4)
        self.assertAlmostEqual(y[0][2], 0., 4)
        self.assertAlmostEqual(y[0][3], 0., 4)
        self.assertAlmostEqual(y[0][4], 0., 4)

        # check other points on the trajectory
        self.assertAlmostEqual(y[1][0], 0.9231163463, 4)
        self.assertAlmostEqual(y[13][0], 0.353454681, 4)
        self.assertAlmostEqual(y[8][1], 0.142837751, 4)
        self.assertAlmostEqual(y[20][1], 0.492844600, 4)
        self.assertAlmostEqual(y[15][2], 0.346376313, 4)
        self.assertAlmostEqual(y[27][2], 0.230837103, 4)
        self.assertAlmostEqual(y[22][3], 0.081296859, 4)
        self.assertAlmostEqual(y[37][3], 0.039501126, 4)
        self.assertAlmostEqual(y[29][4], 0.150075280, 4)
        self.assertAlmostEqual(y[41][4], 0.078591978, 4)
        self.assertAlmostEqual(y[50][0], 0.018315639, 4)
        self.assertAlmostEqual(y[50][1], 0.917958431, 4)
        self.assertAlmostEqual(y[50][2], 0.06372593, 4)
        self.assertAlmostEqual(y[50][3], 0.018207409, 4)
        self.assertAlmostEqual(y[50][4], 0.045518522, 4)




    def test_enforce_non_negativity(self):
        """ Test enforcement of non-negativity during integration """

        # check to make sure that the answer is *incorrect* if we don't enforce
        # nonegativity (ineq_constr=0)
            
        y, t, ypout, t_root, y_root, i_root = daeint(non_neg_res_func, 
                                                     tlist_non_neg,
                                                     y0_non_neg, yp0_non_neg,
                                                     rtol = reltol_non_neg,
                                                     atol = abstol_non_neg,
                                                     ineq_constr=False)


        self.assertAlmostEqual(y[1][0], 0.960000000, 4)
        self.assertAlmostEqual(y[-4][0], -.8800000000, 4)


        # check to make sure that the answer is *correct* if we do enforce
        # nonegativity (ineq_constr=2)

        y, t, ypout, t_root, y_root, i_root = daeint(non_neg_res_func, 
                                                     tlist_non_neg,
                                                     y0_non_neg, yp0_non_neg,
                                                     rtol = reltol_non_neg,
                                                     atol = abstol_non_neg,
                                                     ineq_constr=True)

        self.assertAlmostEqual(y[1][0], 0.960000000, 4)
        self.assertAlmostEqual(y[-4][0], 0.000000, 4)


    def test_redirect_output(self):
        """ Test to make sure we can turn output redirection on and off """
        
        # By default output redirection is off, so we begin by doing an example
        # that should generate errors and making sure that no output is received.
        redir = Utility.Redirector()
        redir.start()
        # This example will generate errors because the maximum number of steps
        # (500) will be passed
        y, t, ypout, t_root, y_root, i_root = daeint(trig_res_func, tlist_trig,
                                                     y0_trig, yp0_trig,
                                                     rtol = reltol_trig,
                                                     atol = abstol_trig,
                                                     max_steps = 7500)
        messages = redir.stop()
        self.assertEqual(len(messages), 0)

        redir = Utility.Redirector()
        redir.start()
        # Now we do the same example again with output redirection off
        y, t, ypout, t_root, y_root, i_root = daeint(trig_res_func, tlist_trig,
                                                     y0_trig, yp0_trig,
                                                     rtol = reltol_trig,
                                                     atol = abstol_trig,
                                                     max_steps = 7500,
                                                     redir_output = False)

        messages = redir.stop()
        self.assertNotEqual(len(messages), 0)


################################################################################

        
suite = unittest.makeSuite(test_daskr)

if __name__ == '__main__':
    unittest.main()
