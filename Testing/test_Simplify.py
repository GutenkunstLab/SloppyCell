from compiler.ast import *
import random
import sets
import unittest

import SloppyCell.ExprManip as ExprManip

x = random.random()
y = random.random()
z = random.random()
f = lambda x: x**4
g = lambda x, y, z: x**4 * y**6 * z**5


class test_Simplify(unittest.TestCase):
    def test_simplify_expr(self):
        cases = ['x', 'x+y', 'x-y', 'x*y', 'x/y', 'x**y', '-x', 'x**-y',
                 'x**(-y + z)', 'f(x)', 'g(x,y,z)', 'x**(y**z)', 
                 '(x**y)**z', 'x**y**z', 'x - (x+y)', '(x+y) - z',
                 'g(x-0+2, y**2 - 0**0, z*y + x/1)', 'x/x', 'x/y',
                 '(x-x)/z', 'x**2 - y/z', 'x+1-1+2-3-x', '0+1*1', 'x-x+y', 
                 '(-2)**2', '-2**2', 'x/y == x/y', 'not True', 'x/x + y/y == 2',
                 '3 + 4 > 6', '3 + (4 > 6)',
                 ]

        for expr in cases: 
            simplified = ExprManip.simplify_expr(expr)
            orig = eval(expr)
            simp = eval(simplified)
            if orig != 0:
                assert abs(orig - simp)/(0.5 * (orig + simp)) < 1e-6
            else:
                assert simp == 0

suite = unittest.makeSuite(test_Simplify)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
