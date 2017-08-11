"""
Functions for simplifying python math expressions.
"""
from compiler.ast import *
import operator
import sets

import AST

# Constants for comparison
_ZERO = Const(0)
_ONE = Const(1)

def simplify_expr(expr):
    """
    Return a simplified version of the expression.
    """
    ast = AST.strip_parse(expr)
    return AST.ast2str(_simplify_ast(ast))

def _simplify_ast(ast):
    """
    Return a simplified ast.

    Current simplifications:
        Special cases for zeros and ones, and combining of constants, in 
            addition, subtraction, multiplication, division.
        Note that at present we only handle constants applied left to right.
          1+1+x -> 2+x, but x+1+1 -> x+1+1.
        x - x = 0
        --x = x
    """
    if isinstance(ast, Name) or isinstance(ast, Const):
        return ast
    elif isinstance(ast, Add) or isinstance(ast, Sub):
        # We collect positive and negative terms and simplify each of them
        pos, neg = [], []
        AST._collect_pos_neg(ast, pos, neg)
        pos = [_simplify_ast(term) for term in pos]
        neg = [_simplify_ast(term) for term in neg]
        
        # We collect and sum the constant values
        values = [term.value for term in pos if isinstance(term, Const)] +\
                [-term.value for term in neg if isinstance(term, Const)] 
        value = sum(values)

        # Remove the constants from our pos and neg lists
        pos = [term for term in pos if not isinstance(term, Const)]
        neg = [term for term in neg if not isinstance(term, Const)]

        new_pos, new_neg = [], []
        for term in pos:
            if isinstance(term, UnarySub):
                new_neg.append(term.expr)
            else:
                new_pos.append(term)
        for term in neg:
            if isinstance(term, UnarySub):
                new_pos.append(term.expr)
            else:
                new_neg.append(term)
        pos, neg = new_pos, new_neg

        # Append the constant value sum to pos or neg
        if value > 0:
            pos.append(Const(value))
        elif value < 0:
            neg.append(Const(abs(value)))

        # Count the number of occurances of each term.
        term_counts = [(term, pos.count(term) - neg.count(term)) for term in
                       pos + neg]
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        term_counts = dict([(str(term), (term, count)) for term, count in 
                            term_counts])

        # We find the first term with non-zero count.
        ii = 0
        for ii, term in enumerate(pos+neg):
            ast_out, count = term_counts[str(term)]
            if count != 0:
                break
        else:
            # We get here if we don't break out of the loop, implying that
            #  all our terms had count of 0
            return _ZERO

        term_counts[str(term)] = (ast_out, 0)
        if abs(count) != 1:
            ast_out = Mul((Const(abs(count)), ast_out))
        if count < 0:
            ast_out = UnarySub(ast_out)

        # And add in all the rest
        for term in (pos+neg)[ii:]:
            term, count = term_counts[str(term)]
            term_counts[str(term)] = (term, 0)
            if abs(count) != 1:
                term = Mul((Const(abs(count)), term))
            if count > 0:
                ast_out = Add((ast_out, term))
            elif count < 0:
                ast_out = Sub((ast_out, term))

        return ast_out
    elif isinstance(ast, Mul) or isinstance(ast, Div):
        # We collect numerator and denominator terms and simplify each of them
        num, denom = [], []
        AST._collect_num_denom(ast, num, denom)
        num = [_simplify_ast(term) for term in num]
        denom = [_simplify_ast(term) for term in denom]
        
        # We collect and sum the constant values
        values = [term.value for term in num if isinstance(term, Const)] +\
                [1./term.value for term in denom if isinstance(term, Const)] 
        # This takes the product of all our values
        value = reduce(operator.mul, values + [1])

        # If our value is 0, the expression is 0
        if not value:
            return _ZERO

        # Remove the constants from our pos and neg lists
        num = [term for term in num if not isinstance(term, Const)]
        denom = [term for term in denom if not isinstance(term, Const)]

        # Here we count all the negative (UnarySub) elements of our expression.
        # We also remove the UnarySubs from their arguments. We'll correct
        #  for it at the end.
        num_neg = 0
        for list_of_terms in [num, denom]:
            for ii, term in enumerate(list_of_terms):
                if isinstance(term, UnarySub):
                    list_of_terms[ii] = term.expr
                    num_neg += 1

        # Append the constant value sum to pos or neg
        if abs(value) != 1:
            num.append(Const(abs(value)))
        if value < 0:
            num_neg += 1

        make_neg = num_neg % 2

        # Count the number of occurances of each term.
        term_counts = [(term, num.count(term) - denom.count(term)) for term in
                       num + denom]
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        term_counts = dict([(str(term), (term, count)) for term, count in 
                            term_counts])

        nums, denoms = [], []
        # We walk through terms in num+denom in order, so we rearrange a little
        #  as possible.
        for term in num+denom:
            term, count = term_counts[str(term)]
            # Once a term has been done, we set its term_counts to 0, so it
            #  doesn't get done again.
            term_counts[str(term)] = (term, 0)
            if abs(count) > 1:
                term = Power((term, Const(abs(count))))
            if count > 0:
                nums.append(term)
            elif count < 0:
                denoms.append(term)

        # We return the product of the numerator terms over the product of the
        #  denominator terms
        out = AST._make_product(nums)
        if denoms:
            denom = AST._make_product(denoms)
            out = Div((out, denom))

        if make_neg:
            out = UnarySub(out)

        return out
    elif isinstance(ast, Power):
        # These cases all have a left and a right, so we group them just to
        #  avoid some code duplication.
        power = _simplify_ast(ast.right)
        base = _simplify_ast(ast.left)

        if power == _ZERO:
            # Anything, including 0, to the 0th power is 1, so this
            #  test should come first
            return _ONE
        if base == _ZERO or base == _ONE or power == _ONE:
            return base
        elif isinstance(base, Const) and\
                isinstance(power, Const):
            return Const(base.value**power.value)
        # Getting here implies that no simplifications are possible, so just
        #  return with simplified arguments
        return Power((base, power))
    elif isinstance(ast, UnarySub):
        simple_expr = _simplify_ast(ast.expr)
        if isinstance(simple_expr, UnarySub):
            # Case --x
            return _simplify_ast(simple_expr.expr)
        elif isinstance(simple_expr, Const):
            if simple_expr.value == 0:
                return Const(0)
            else:
                return Const(-simple_expr.value)
        else:
            return UnarySub(simple_expr)
    elif isinstance(ast, UnaryAdd):
        simple_expr = _simplify_ast(ast.expr)
        return simple_expr
    elif isinstance(ast, list):
        simple_list = [_simplify_ast(elem) for elem in ast]
        return simple_list
    elif isinstance(ast, tuple):
        return tuple(_simplify_ast(list(ast)))
    elif AST._node_attrs.has_key(ast.__class__):
        # Handle node types with no special cases.
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            if isinstance(attr, list):
                for ii, elem in enumerate(attr):
                    attr[ii] = _simplify_ast(elem)
            else:
                setattr(ast, attr_name, _simplify_ast(attr))
        return ast
    else:
        return ast
