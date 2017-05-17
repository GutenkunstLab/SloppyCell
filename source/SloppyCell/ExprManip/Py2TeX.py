from compiler.ast import *
import os

import AST

# This is just and istance of Mul to use when we group numerators and 
#  denominators
_EMPTY_MUL = Mul((None, None))

def dict2TeX(d, name_dict, lhs_form='%s', split_terms=False, simpleTeX=False):
    lines = []
    for lhs, rhs in d.items():
        if split_terms:
            ast = AST.strip_parse(rhs)
            pos, neg = [], []
            AST._collect_pos_neg(ast, pos, neg)
            try:
                lhsTeX = lhs_form % expr2TeX(lhs, name_dict=name_dict)
            except TypeError:
                lhsTeX = lhs_form
            rhsTeX = _ast2TeX(pos[0], name_dict=name_dict)
            lines.append(r'$ %s $ &=& $ %s $\\' % (lhsTeX, rhsTeX))
                         
            for term in pos[1:]:
                TeXed = _ast2TeX(term, name_dict=name_dict)
                lines.append(r' & & $ + \, %s $\\' % TeXed)
            for term in neg:
                TeXed = _ast2TeX(term, name_dict=name_dict)
                lines.append(r' & & $ - \, %s $\\' % TeXed)
        else:
            lhsTeX = lhs_form % expr2TeX(lhs, name_dict=name_dict)
            rhsTeX = expr2TeX(rhs, name_dict=name_dict)
            lines.append(r'$ %s $ & = & $ %s $\\' % (lhsTeX, rhsTeX))

        if not simpleTeX:
            # Force a space between TeX'd entries
            lines[-1] = '%s[5mm]' % lines[-1]

    all = os.linesep.join(lines)

    if not simpleTeX:
        all = all.replace(r'\frac{', r'\tabfrac{')
        # This makes the fractions look much nicer in the tabular output. See
        #  http://www.texnik.de/table/table.phtml#fractions
        lines = [r'\providecommand{\tabfrac}[2]{%',
                 r'   \setlength{\fboxrule}{0pt}%',
                 r'   \fbox{$\frac{#1}{#2}$}}',
                 r'\begin{longtable}{lll}'] + [all] + [r'\end{longtable}']
        all = os.linesep.join(lines)

    return all

def expr2TeX(expr, name_dict={}):
    """
    Return a TeX version of a python math expression.

    name_dict: A dictionary mapping variable names used in the expression to
        preferred TeX expressions.
    """
    ast = AST.strip_parse(expr)
    return _ast2TeX(ast, name_dict=name_dict)

def _ast2TeX(ast, outer=AST._FARTHEST_OUT, name_dict={}, 
             adjust=0):
    """
    Return a TeX version of an AST.

    outer: The AST's 'parent' node, used to determine whether or not to 
        enclose the result in parentheses. The default of _FARTHEST_OUT will
        never enclose the result in parentheses.

    name_dict: A dictionary mapping variable names used in the expression to
        preferred TeX expressions.

    adjust: A numerical value to adjust the priority of this ast for
        particular cases. For example, the denominator of a '/' needs 
        parentheses in more cases than does the numerator.
    """
    if isinstance(ast, Name):
        # Try to get a value from the name_dict, defaulting to ast.name if
        #  ast.name isn't in name_dict
        out = name_dict.get(ast.name, ast.name)
    elif isinstance(ast, Const):
        out = str(ast.value)
    elif isinstance(ast, Add):
        out = '%s + %s' % (_ast2TeX(ast.left, ast, name_dict),
                           _ast2TeX(ast.right, ast, name_dict))
    elif isinstance(ast, Sub):
        out = '%s - %s' % (_ast2TeX(ast.left, ast, name_dict),
                           _ast2TeX(ast.right, ast, name_dict, adjust = 1))
    elif isinstance(ast, Mul) or isinstance(ast, Div):
        # We collect all terms numerator and denominator
        nums, denoms = [], []
        AST._collect_num_denom(ast, nums, denoms)
        # _EMPTY_MUL ensures that parentheses are done properly, since every
        #  element is now the child of a Mul
        lam_func = lambda arg: _ast2TeX(arg, _EMPTY_MUL, name_dict)
        nums = [lam_func(term) for term in nums]
        if denoms:
            denoms = [lam_func(term) for term in denoms]
            out = r'\frac{%s}{%s}' % (r' \cdot '.join(nums),
                                      r' \cdot '.join(denoms))
        else:
            out = r' \cdot '.join(nums)
    elif isinstance(ast, Power):
        out = '{%s}^{%s}' % (_ast2TeX(ast.left, ast, name_dict, adjust = 1),
                             _ast2TeX(ast.right, ast, name_dict))
    elif isinstance(ast, UnarySub):
        out = '-%s' % _ast2TeX(ast.expr, ast, name_dict)
    elif isinstance(ast, UnaryAdd):
        out = '+%s' % _ast2TeX(ast.expr, ast, name_dict)
    elif isinstance(ast, CallFunc):
        lam_func = lambda arg: _ast2TeX(arg, name_dict=name_dict)
        name = lam_func(ast.node)
        args = [lam_func(arg) for arg in ast.args]
        if name == 'sqrt' and len(args) == 1:
            # Special case
            out = r'\sqrt{%s}' % args[0]
        else:
            out = r'\operatorname{%s}\left(%s\right)' % (name, r',\,'.join(args))
    elif isinstance(ast, Or) or isinstance(ast, And):
        out = r'\operatorname{%s}' % (str(ast))
        
    if AST._need_parens(outer, ast, adjust):
        return out
    else:
        return r'\left(%s\right)' % out
