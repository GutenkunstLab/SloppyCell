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

        # Append the constant value sum to pos or neg
        if value > 0:
            pos.append(Const(value))
        elif value < 0:
            neg.append(Const(-value))

        # Count the number of occurances of each term.
        term_counts = [(term, pos.count(term) - neg.count(term)) for term in
                       pos + neg]
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        term_counts = dict([(str(term), (term, count)) for term, count in 
                            term_counts if count != 0])

        # If all our terms have canceled
        if not term_counts:
            return _ZERO

        # We deal with the first term
        temp, (ast_out, count) = term_counts.popitem()
        if count == -1:
            ast_out = UnarySub(ast_out)
        elif count != 1:
            ast_out = Mul((Const(count), ast_out))

        # And add in all the rest
        for term, count in term_counts.values():
            if count > 1:
                term = Mul((Const(count), term))
            elif count < -1:
                term = Mul((Const(-count), term))
            if count > 0:
                ast_out = Add((ast_out, term))
            else:
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
        if values:
            value = reduce(operator.mul, values)
        else:
            value = 1

        # If our value is 0, the expression is 0
        if not value:
            return _ZERO

        # Remove the constants from our pos and neg lists
        num = [term for term in num if not isinstance(term, Const)]
        denom = [term for term in denom if not isinstance(term, Const)]

        # Append the constant value sum to pos or neg
        if abs(value) != 1:
            num.append(Const(abs(value)))
        if value < 0:
            make_neg = True
        else:
            make_neg = False


        # Count the number of occurances of each term.
        term_counts = [(term, num.count(term) - denom.count(term)) for term in
                       num + denom]
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        term_counts = dict([(str(term), (term, count)) for term, count in 
                            term_counts if count != 0])

        nums, denoms = [], []
        for term, count in term_counts.values():
            if abs(count) > 1:
                term = Power((term, Const(abs(count))))
            if count > 0:
                nums.append(term)
            else:
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
    else:
        # Handle node types with no special cases.
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            if isinstance(attr, list):
                for ii, elem in enumerate(attr):
                    attr[ii] = _simplify_ast(elem)
            else:
                setattr(ast, attr_name, _simplify_ast(attr))
        return ast
