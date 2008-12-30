import random
import sets
import unittest

import SloppyCell.ExprManip as ExprManip

x = random.random()
y = random.random()
z = random.random()
p = [random.random(), random.random()]
f = lambda x: x**4
g = lambda x, y, z: x**4 * y**6 * z**5

class test_Substitution(unittest.TestCase):
    def test_sub_for_var(self):
        cases = [('x', 'x', 'y', 
                  'y'),
                 ('x + 1', 'x', 'y', 
                  'y + 1'),
                 ('f(x)**2 - y + z', 'x', 'y*z/(x-2)', 
                  'f(y*z/(x-2))**2 - y + z'),
                 ('x**2', 'x', 'y+2', 
                  '(y+2)**2'),
                 ('p[0]**x', 'x', 'y+2', 
                  'p[0]**(y+2)'),
                 ('g(x, y, f(z)) != y+2 and z == y', 'z', 'y',
                  '(g(x, y, f(y)) != (y + 2)) and (y == y)'),
                 ]

        for expr, out_var, in_expr, answer in cases:
            subbed = ExprManip.sub_for_var(expr, out_var, in_expr)
            assert eval(answer) == eval(subbed)

    def test_sub_for_vars(self):
        cases = [('x', {'x':'y'},
                  'y'),
                 ('x + 1', {'x':'y'}, 
                  'y + 1'),
                 ('f(x)**2 - y + z', {'x':'y*z/(x-2)'}, 
                  'f(y*z/(x-2))**2 - y + z'),
                 ('x**2', {'x':'y+2'}, 
                  '(y+2)**2'),
                 ('p[0]**x', {'x':'y+2'}, 
                  'p[0]**(y+2)'),
                # In these cases substituting one-by-one will fail
                 ('x*y', {'x':'y', 'y':'z'},
                  'y*z'),
                 ('z*y', {'z':'y', 'y':'z'},
                  'y*z'),
                 ]

        for expr, mapping, answer in cases:
            subbed = ExprManip.sub_for_vars(expr, mapping)
            assert eval(answer) == eval(subbed)

    def test_sub_for_comps(self):
        cases = [('not x < 3 and True', {'x < 3': 'True'}, False),
                 ('not x < 3 and x < 3', {'x < 3': 'True'}, False),
                 ('x < 3 and x < 3', {'x < 3': 'True'}, True),
                 ('x < 3 and y > 4', {'x < 3': 'True', 'y > 4':'False'}, False),
                 ('x < 3 and not y > 4', {'x < 3': 'True', 'y > 4':'False'}, 
                  True),]
        for expr, mapping, result in cases:
            subbed = ExprManip.sub_for_comps(expr, mapping)
            assert eval(subbed) == result


    def test_sub_for_func(self):
        cases = [('f(x)', 'f', 'y', 'y+1',
                  'x+1'),
                 ('f(x)**2', 'f', 'y', 'y+1',
                  '(x+1)**2'),
                 ('f(g(x, y, z))**2', 'f', 'y', 'y+1',
                  '(g(x,y,z)+1)**2'),
                 ('g(f(x), y, z)**2', 'f', 'y', 'y+1',
                  'g(x+1,y,z)**2'),
                 ('g(f(x), p[0], z)**2', 'f', 'y', 'y+1',
                  'g(x+1,p[0],z)**2'),
                 ('g(z, y, f(x)) == y+2 and f(x) + f(y) != f(x)', 
                  'g', ('x','y', 'z'), 'f(x + y/z)',
                  '(f(z + y/f(x)) == (y + 2)) and (f(x) + f(y) != f(x))'),
                 ]

        for expr, func_name, func_vars, func_expr, answer in cases:
            subbed = ExprManip.sub_for_func(expr, func_name, func_vars,
                                               func_expr)
            assert eval(answer) == eval(subbed)

    def test_var_args_subs(self):
        subbed = ExprManip.sub_for_func('or_func(x,y,z)', 'or_func', '*',
                                        'a or b')
        x,y,z = True,False,False
        assert eval(subbed)
        x,y,z = False,False,False
        assert not eval(subbed)

        self.assertRaises(ValueError, ExprManip.sub_for_func,
                          'or_func(x,y,z)', 'or_func', '*', 'a + b')

suite = unittest.makeSuite(test_Substitution)

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite)
