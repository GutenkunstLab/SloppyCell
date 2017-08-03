import random
import sets
import unittest

from math import *

import SloppyCell.ExprManip as ExprManip

# y is larger up to avoid negative values under sqrts in tests.
x, y, z = 0.71, 8.9, 0.02
f = lambda x: x**4
f_0 = lambda x: 4*x**3
g = lambda x, y, z: x**4 * y**6 * z**5
g_0 = lambda x, y, z: 4*x**3 * y**6 * z**5
g_1 = lambda x, y, z: x**4 * 6*y**5 * z**5
g_2 = lambda x, y, z: x**4 * y**6 * 5*z**4


class test_Differentiation(unittest.TestCase):
    def test_diff_expr(self):
        cases = [('x', 'x'),
                 ('+x', 'x'),
                 ('y', 'x'),
                 ('x+y', 'x'),
                 ('x*y', 'x'),
                 ('x/y', 'x'),
                 ('x**y', 'x'),
                 ('x**x', 'x'),
                 ('y**x', 'x'),
                 ('acos(x)', 'x'),
                 ('asin(x)', 'x'),
                 ('atan(x)', 'x'),
                 ('cos(x)', 'x'),
                 ('cosh(x)', 'x'),
                 ('exp(x)', 'x'),
                 ('log(x)', 'x'),
                 ('log10(x)', 'x'),
                 ('sin(x)', 'x'),
                 ('sinh(x)', 'x'),
                 ('sqrt(x)', 'x'),
                 ('tan(x)', 'x'),
                 ('tanh(x)', 'x'),
                 ('pow(x, 2)', 'x'),
                 ('pow(y, x)', 'x'),
                 ('y**2 - 4*x*z', 'x'),
                 ('(-y + sqrt(y**2 - 4*x*z))/(2*x)', 'x'),
                 ('(-y + sqrt(y**2 - 4*x*z))/(2*x)', 'y'),
                 ('(-y + sqrt(y**2 - 4*x*z))/(2*x)', 'z'),
                 ('sqrt(y - x)', 'x'),
                 ('f(x)' ,'x'),
                 ('f(x - y*x + cos(x*z))' ,'x'),
                 ('g(x - y*x + cos(x*z), y, x)' ,'x'),
                 ('g(x - y*x + cos(x*z), y, x)' ,'y'),
                 ('g(x - y*x + cos(x*z), y**z, x - 14*z)' ,'y'),
                 ('g(x - y*x + cos(x*z), y, x)' ,'y'),
                 ('g(x*y, x*y**2, y)', 'z')
                 ]

        for expr, wrt in cases: 
            d = ExprManip.diff_expr(expr, wrt)
            ad = eval(d)
            fd = self._num_diff(expr, wrt, x=x, y=y, z=z)
            # We test that our numeric and analytic derivatives differ by less
            #  than 0.1%
            if ad != 0:
                assert abs(fd - ad)/(0.5*(ad + fd)) < 1e-3
            else:
                assert ad == fd

    def _num_diff(self, expr, wrt, delta = 1e-3, **args):
        """
        Take the numerical derivative of expr. Variables used in expr should be
        passed in as keyword arguments.
        """
        # We do our evals using args as the locals dictionary, so we get
        #  whatever values are passed in.
        wrt_val = args[wrt]
        args[wrt] = wrt_val*(1-delta/2.0)
        fmin = eval(expr, globals(), args)
        args[wrt] = wrt_val*(1+delta/2.0)
        fplus = eval(expr, globals(), args)
    
        return (fplus - fmin)/(wrt_val * delta)

suite = unittest.makeSuite(test_Differentiation)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
