import os
import unittest

import SloppyCell.ExprManip as ExprManip

lines = []

class test_Py2TeX(unittest.TestCase):
    def test_expr2Tex(self):
        lines.append(r'\section{Basic expressions}')
        lines.append(r'\begin{longtable}{rl}')
        cases = ['x', 'x+y', 'x-y', 'x*y', 'x/y', 'x**y', '-x', 'x**-y',
                 'x**(-y + z)', 'f(x)', 'g(x,y,z)', 'x**(y**z)', 
                 '(x**y)**z', 'x**y**z', 'x - (x+y)', '(x+y) - z',
                 'g(x-0+2, y**2 - 0**0, z*y + x/1)', 'sqrt(x+y-sqrt(z/x))',
                 ]

        for expr in cases: 
            TeXed = ExprManip.expr2TeX(expr)
            line = r'{\tt %s} & $ %s $ \\' % (expr, TeXed)
            lines.append(line)

        lines.append(r'\end{longtable}')

    def test_name_dict(self):
        name_dict = {'alpha': r'\alpha',
                   'beta': r'\beta',
                   'betarho': r'\beta\rho',
                   'chi': r'\lambda',
                   'tau': r'\tau'
                   }
        cases = ['alpha',
                 'alpha**beta',
                 'chi(betarho, tau)',
                 'sqrt(4+sqrt(alpha))',
                 ]

        for expr in cases: 
            TeXed = ExprManip.expr2TeX(expr, name_dict)

suite = unittest.makeSuite(test_Py2TeX)

if __name__ == '__main__':
    lines.append(r'\documentclass{article}')
    lines.append(r'\usepackage{amsmath}')
    lines.append(r'\usepackage{fullpage}')
    lines.append(r'\usepackage{longtable}')
    lines.append(r'\begin{document}')
    unittest.TextTestRunner(verbosity=2).run(suite)
    lines.append(r'\end{document}')
    f = file('test_Py2TeX.tex', 'w')
    f.write(os.linesep.join(lines))
    f.close()
