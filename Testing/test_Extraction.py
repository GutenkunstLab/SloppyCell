import sets
import unittest

import SloppyCell.ExprManip as ExprManip

class test_Extraction(unittest.TestCase):
    def test_extract_vars(self):
        cases = [('x', ['x']),
                 ('x + y', ['x', 'y']),
                 ('x * y', ['x', 'y']),
                 ('x / y', ['x', 'y']),
                 ('x**2', ['x']),
                 ('x**y', ['x', 'y']),
                 ('f(x)', ['x']),
                 ('f(x, y)', ['x', 'y']),
                 ('f(x + z/x, y)', ['x', 'y', 'z']),
                 ('x**2 + f(y)', ['x', 'y']),
                 ('x + y**2 + z**(a + b) + z**f.g.h(a + b + sqrt(c))', 
                  ['x', 'y', 'z', 'a', 'b', 'c']),
                 ('x < y', ['x', 'y']),
                 ('(x < y) and (x == 2)', ['x', 'y']),
                 ]

        for expr, vars in cases:
            assert ExprManip.extract_vars(expr) == sets.Set(vars)

    def test_extract_funcs(self):
        cases = [('g(x)', [('g', 1)]),
                 ('f(g(x))', [('f', 1), ('g', 1)]),
                 ('f(g(x)/2, a, b(x, y))', 
                  [('f', 3), ('g', 1), ('b', 2)]),
                 ('a + g(x)', [('g', 1)]),
                 ('a + g(x) + (a + f(x))', [('f', 1), ('g', 1)]),
                 ('f(g(1 + h(a))/2) + j(q**k(x)) + q.r[1] + s[1]', 
                  [('f', 1), ('g', 1), ('h', 1), ('j', 1), ('k', 1)]),
                 ('(g(x) + 2)**f(y)', [('f', 1), ('g', 1)]),
                 ('g(x)**f(y)', [('f', 1), ('g', 1)]),
                 ('g(x) > f(y)', [('f', 1), ('g', 1)]),
                 ('g(x) and f(y)', [('f', 1), ('g', 1)]),
                 ('g(x) and not f(y)', [('f', 1), ('g', 1)])
                 ]

        for expr, funcs in cases:
            assert ExprManip.extract_funcs(expr) == sets.Set(funcs)

    def test_extract_comps(self):
        cases = [('x == 3', ['x == 3']),
                 ('x == 3 and y == 4', ['x == 3', 'y == 4']),
                 ('x < 3 and not y > 4', ['x < 3', 'y > 4']),
                 ('x < 3 + (y > 4)', ['x < (3 + (y > 4))', 'y > 4']),
                 ]

        for expr, comps in cases:
            result = ExprManip.extract_comps(expr) 
            assert result == sets.Set(comps)

suite = unittest.makeSuite(test_Extraction)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
