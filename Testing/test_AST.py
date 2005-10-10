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
