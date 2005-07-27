import os
import unittest

import scipy

from SloppyCell.ReactionNetworks import *

_TEX_DIR = 'TeX'

class test_TeX(unittest.TestCase):
    def setUp(self):
        if not os.path.isdir(_TEX_DIR):
            os.mkdir(_TEX_DIR)

    def test_basic_TeX(self):
        """Test generation of individual TeX strings"""
        snippets = ['x',
                    'x + y', 
                    'x * y',
                    'x / y',
                    'x**2',
                    'x**y',
                    'f(x)',
                    'f(x, y)',
                    'f(x + z/x.y, y)',
                    'x**2 + f(y)',
                    'z**(a + b) + z**f.g.h(a + b + sqrt(c))',
                    'g(a)**2',
                    'f.g(x)',
                    'f(g(x))',
                    'f(g(x)/2, a, b(x))',
                    'a + g(x)',
                    'a + g(x) + (a + f(x))',
                    'f(g(1 + h.i(a))/2) + j(q**k(x))',
                    '(g(x) + 2)**f(y)',
                    'g(x)**f(y)',
                    'f(b, a) + a(q) + a.b(x) + a**x + y**a',
                    'f.b(x + a.b)',
                    'g(a + b, c)',
                    'f(g(a + b))',
                    'f.g.h(c)',
                    'f.g(c)',
                    'f.g(c)**2',
                    'f(sqrt(q))**3',
                    '4**2**3',
                    'alpha',
                    'alpha**beta',
                    'chi(betarho, tau)',
                    'sqrt(4+sqrt(alpha))',
                    ]

        mapping = {'alpha': r'\alpha',
                   'beta': r'\beta',
                   'betarho': r'\beta\rho',
                   'chi': r'\lambda',
                   'tau': r'\tau'
                   }

        lines = []
        lines.append(r'\documentclass{article}')
        lines.append(r'\usepackage{amsmath}')
        lines.append(r'\begin{document}')
        lines.append(r'\begin{tabular}{r|l}')
        lines.append(r'Python code & TeX form\\')
        lines.append(r'\hline')

        for snip in snippets:
            tex = IO.expr_to_TeX(snip, mapping)
            lines.append('{\\tt %s} & $%s$ \\\\' % (snip, tex))
            lines.append('')

        lines.append(r'\end{tabular}')
        lines.append(r'\end{document}')
        f = file(os.path.join(_TEX_DIR, 'basic.tex'), 'w')
        f.write(os.linesep.join(lines))
        f.close()

suite = unittest.makeSuite(test_TeX, 'test')
message = 'Results of TeX output testing are in the directory "%s".' % _TEX_DIR

if __name__ == '__main__':
    unittest.main()
