import os

import symbolic

def net_eqns_to_TeX_file(net, filename):
    lines = []
    lines.append(r'\documentclass{article}')
    lines.append(r'\usepackage{amsmath}')
    lines.append(r'\usepackage{fullpage}')
    lines.append(r'\usepackage{longtable}')
    lines.append(r'\begin{document}')
    lines.append(_net_eqns_to_TeX(net))
    lines.append(r'\end{document}')

    f = file(filename, 'w')
    f.write(os.linesep.join(lines))
    f.close()

def _net_eqns_to_TeX(net):
    """
    Return a string that contains the longtable-bound TeX'd equations for the network
    """
    name_dict = dict([(id, '\\mathrm{%s}' % net.get_var_name(id))
                      for id in net.variables.keys()])
    species_dict = dict([(id, '\\left[\\mathrm{%s}\\right]'
                          % net.get_var_name(id))
                         for id in net.species.keys()])
    name_dict.update(species_dict)
    lines = []
    # This makes the fractions look much nicer in the tabular output. See
    #  http://www.texnik.de/table/table.phtml#fractions
    lines.append(r'\newcommand{\tabfrac}[2]{%')
    lines.append(r'   \setlength{\fboxrule}{0pt}%')
    lines.append(r'   \fbox{$\frac{#1}{#2}$}}')
    lines.append(r'\begin{longtable}{rcl}')
    for id, rhs in net.diffEqRHS.items():
        texRHS = expr_to_TeX(rhs, name_dict, longtable=True)
        line = '$ \\frac{d\\,%s}{dt}$ &=& %s' % (name_dict.get(id), texRHS)
        # Use the tabfrac and force a space between the equations
        lines.append(line.replace(r'\frac{', r'\tabfrac{') + '[5mm]')
        lines.append('')
    for id, rhs in net.assignmentRules.items():
        texRHS = expr_to_TeX(rhs, name_dict, longtable=True)
        line = '$ %s $ &=& %s' % (name_dict.get(id), texRHS)
        lines.append(line.replace(r'\frac{', r'\tabfrac{') + '[5mm]')
        lines.append('')

    lines.append(r'\end{longtable}')

    return os.linesep.join(lines)


def expr_to_TeX(input, name_dict={}, longtable=False):
    """
    Output the TeX representation of a python expression.

    name_dict: an optional dictionary mapping names in the python expression
               to TeX'd names
    longtable: if True, the top level of an expression is broken into
                 a set of terms suitable for inclusion into a table
    """
    ast = symbolic.string2ast(input)
    # Note that we drop parentheses if they exist around the entire expression
    if ast[0] == 'arith_expr' and  longtable:
        lines = []
        lines.append('$ %s $\\\\' % _ast_to_TeX(_drop_parens(ast[1]),
                                                name_dict))
        for op_index in range(2, len(ast), 2):
            dropped = _drop_parens(ast[op_index+1])
            lines.append(' & & $ %s\\,%s $ \\\\' % (ast[op_index][1], 
                                               _ast_to_TeX(dropped, name_dict)))
        lines[-1] = '%s' % lines[-1]
        return os.linesep.join(lines)
    elif longtable:
        return ' $ %s $ \\\\' % (_ast_to_TeX(_drop_parens(ast), name_dict))
    else:
        return _ast_to_TeX(_drop_parens(ast), name_dict)

def _drop_parens(ast):
    # Drops the parentheses from an AST, if they exist
    if ast[0] == 'atom':
        return ast[2]
    else:
        return ast


def _ast_to_TeX(term, name_dict = {}):
    """
    Return the TeX version of an ast tree.

    If name_dict is not empty, variable names are substituded by their 
    corresponding values in the dictionary.
    """
    if term[0] == 'NAME':
        # Try to get a name from the name_dict, but default to just term[1]
        return name_dict.get(term[1], term[1])
    elif term[0] == 'NUMBER':
        return term[1]
    elif term[0] == 'factor':
        if term[1][0] == 'PLUS':
            # We drop extraneous plus signs
            return _ast_to_TeX(term[2], name_dict)
        elif term[1][0] == 'MINUS':
            return '-%s' % _ast_to_TeX(term[2], name_dict)
    elif term[0] == 'arith_expr':
        out = _ast_to_TeX(term[1], name_dict)
        for op_index in range(2, len(term), 2):
            out = '%s %s %s' % (out, term[op_index][1], 
                                _ast_to_TeX(term[op_index+1], name_dict))
        return out
    elif term[0] == 'term':
        # We collect all the numerators and denominators together, so x/y*a/b
        #  comes out as \frac{x*a}{y*b}
        nums = [term[1]]
        denoms = []
        for op_index in range(2, len(term), 2):
            if term[op_index][0] == 'STAR':
                nums.append(term[op_index + 1])
            elif term[op_index][0] == 'SLASH':
                denoms.append(term[op_index+1])
            else:
                print 'Unexpected operator in term %s!' % term

            numerator = r' \times '.join([_ast_to_TeX(t, name_dict)
                                          for t in nums])
            if len(denoms) == 0:
                out = numerator
            elif len(denoms) == 1:
                if denoms[0][0] == 'atom':
                    denominator = _ast_to_TeX(denoms[0][2], name_dict)
                else:
                    denominator = _ast_to_TeX(denoms[0], name_dict)
                out = '\\frac{%s}{%s}' % (numerator, denominator)
            else:
                denominator = r' \times '.join([ast_to_TeX(t, name_dict)
                                                for t in denoms])
                out = '\\frac{%s}{%s}' % (numerator, denominator)
        return out
    elif term[0] == 'atom':
        return '\\left( %s \\right)' % _ast_to_TeX(term[2], name_dict)
    elif term[0] == 'arglist' or term[0] == 'subscriptlist':
        tex_args = [_ast_to_TeX(arg, name_dict) for arg in term[1::2]]
        return r',\, '.join(tex_args)
    # Now we treat the 'power' cases. These are exponents and function calls and
    #  variables of the form a.b.c  The parser doesn't give really nice results
    #  for these cases, so we have to be tricky. Note that the order of the
    #  cases here does matter.
    #
    # These are exponents
    elif term[0] == 'power' and term[-2][0] == 'DOUBLESTAR':
        out = '%s^{%s}' % (_ast_to_TeX(term[:-2], name_dict), 
                           _ast_to_TeX(term[-1], name_dict))
        return out
    # These are function calls
    elif term[0] == 'power' and term[-1][0] == 'trailer' and\
            term[-1][1][0] == 'LPAR':

        # We do sqrt specially
        if term[1][1] == 'sqrt':
            return '\\sqrt{%s}' % _ast_to_TeX(term[-1][2], name_dict)

        out = '\\operatorname{%s}\\left(%s\\right)' % \
                (_ast_to_TeX(term[:-1], name_dict), 
                 _ast_to_TeX(term[-1][2], name_dict))
        return out
    # This is a[1, 3, a]
    # (We don't handle slicing right now)
    elif term[0] == 'power' and term[-1][0] == 'trailer' and\
            term[-1][1][0] == 'LSQB':
        out = '%s\\left[%s\\right]' % \
                (_ast_to_TeX(term[:-1], name_dict), 
                 _ast_to_TeX(term[-1][2], name_dict))
        return out
    # This is a.b.c
    elif term[0] == 'power' and term[-1][0] == 'trailer':
        return '.'.join([term[1][1]] + [t[2][1] for t in term[2:]])
    # This is a kludge to handle pass-downs from the other cases
    elif term[0] == 'power':
        return _ast_to_TeX(term[-1], name_dict)
