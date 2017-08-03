from compiler.ast import *
import random
import sets
import unittest

from SloppyCell.ExprManip import strip_parse, ast2str
import SloppyCell.ExprManip.AST as AST

x = random.random()
y = random.random()
z = random.random()
f = lambda x: x**4
g = lambda x, y, z: x**4 * y**6 * z**5


class test_AST(unittest.TestCase):
    def test_ast2str(self):
        cases = ['x', 'x+y', 'x-y', 'x*y', 'x/y', 'x**y', '-x', 'x**-y',
                 'x**(-y + z)', 'f(x)', 'g(x,y,z)', 'x**(y**z)', 
                 '(x**y)**z', 'x**y**z', 'x - (x+y)', '(x+y) - z',
                 'g(x-0+2, y**2 - 0**0, z*y + x/1)', 'x/x', 'x/y',
                 '(x-x)/z', 'x**2 - y/z', 'x+1-1+2-3-x', '0+1*1', 'x-x+y',
                 '(-2)**2', '-2**2', 'x < 0.5', 'x + y < 0.8', 'x > 0.5',
                 'x+y < z - x', '(f(x) < 0.5) == False', 'x == x', 
                 'x + y == 2*x - x + y', 'True and False', 
                 'True and False or True', 'True or False', 
                 '(True and not False) or True', 'not (True and False)',
                 'x == x and y == y', 'x - x == 0 or y - x != 0']
        for expr in cases: 
            run = ast2str(strip_parse(expr))
            orig = eval(expr)
            out = eval(run)
            if orig != 0:
                assert abs(orig - out)/(0.5 * (orig + out)) < 1e-6
            else:
                assert out == 0

    def test__collect_num_denom(self):
        cases = [(strip_parse('1'), (['1'], [])),
                 (strip_parse('1/2'), (['1'], ['2'])),
                 (strip_parse('1/2*3'), (['1', '3'], ['2'])),
                 (strip_parse('1/(2*3)'), (['1'], ['2', '3'])),
                 (strip_parse('1/(2/3)'), (['1', '3'], ['2'])),
                 (Mul((Mul((Const(1), Const(2))), Mul((Const(3), Const(4))))),
                  (['1', '2', '3', '4'], [])),
                 (Mul((Mul((Const(1), Const(2))), Div((Const(3), Const(4))))),
                  (['1', '2', '3'], ['4'])),
                 (Mul((Div((Const(1), Const(2))), Div((Const(3), Const(4))))),
                  (['1', '3'], ['2', '4'])),
                 (Div((Div((Const(1), Const(2))), Div((Const(3), Const(4))))),
                  (['1', '4'], ['2', '3'])),
                 ]
        for ast, (nums, denoms) in cases: 
            n, d = [], []
            AST._collect_num_denom(ast, n, d)
            n = [ast2str(term) for term in n]
            d = [ast2str(term) for term in d]
            assert sets.Set(nums) == sets.Set(n)
            assert sets.Set(denoms) == sets.Set(d)

    def test__collect_pos_neg(self):
        cases = [(strip_parse('1'), (['1'], [])),
                 (strip_parse('1-2'), (['1'], ['2'])),
                 (strip_parse('1-2+3'), (['1', '3'], ['2'])),
                 (strip_parse('1-(2+3)'), (['1'], ['2', '3'])),
                 (strip_parse('1-(2-3)'), (['1', '3'], ['2'])),
                 (strip_parse('(1-2)-(3-4)'), (['1', '4'], ['2', '3'])),
                 (Add((Add((Const(1), Const(2))), Add((Const(3), Const(4))))),
                  (['1', '2', '3', '4'], [])),
                 (Add((Add((Const(1), Const(2))), Sub((Const(3), Const(4))))),
                  (['1', '2', '3'], ['4'])),
                 (Add((Sub((Const(1), Const(2))), Sub((Const(3), Const(4))))),
                  (['1', '3'], ['2', '4'])),
                 (Sub((Sub((Const(1), Const(2))), Sub((Const(3), Const(4))))),
                  (['1', '4'], ['2', '3'])),
                 ]
        for ast, (poss, negs) in cases: 
            p, n = [], []
            AST._collect_pos_neg(ast, p, n)
            p = [ast2str(term) for term in p]
            n = [ast2str(term) for term in n]
            assert sets.Set(poss) == sets.Set(p)
            assert sets.Set(negs) == sets.Set(n)

suite = unittest.makeSuite(test_AST)

if __name__ == '__main__':
    unittest.main()
