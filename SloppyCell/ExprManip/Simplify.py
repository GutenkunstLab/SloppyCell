"""
Functions for simplifying python math expressions.
"""
from ast import *
import operator
from SloppyCell.ExprManip  import AST
import functools

# Constants for comparison
_ZERO = Constant(value=0)
_ONE = Constant(value=1)


def simplify_expr(expr):
    """
    Return a simplified version of the expression.
    """
    tree = AST.strip_parse(expr)
    simplify_ast = _simplify_ast(tree)
    # print(simplify_ast)
    return AST.ast2str(simplify_ast)

def get_count_from_ast(ast_list, term):
    count = 0
    for item in ast_list:
        if AST.ast2str(item)==AST.ast2str(term):
            count += 1
    return count  
    
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
    # print(dump(ast))
    # try:
    #     # print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",dump(ast))
    #     # print("eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",unparse(ast))
    # except Exception as e:
    #     # print("exception", e)
    if isinstance(ast, Name) or isinstance(ast, Constant):
        return ast
    elif isinstance(ast, BinOp) and (isinstance(ast.op, Add) or isinstance(ast.op, Sub)):
        # print("addition +++++++++++++++++++++++++++++++++++")
        # We collect positive and negative terms and simplify each of them
        pos, neg = [], []
        AST._collect_pos_neg(ast, pos, neg)
        # print("first ---------------")
        # print(pos, neg)
        pos = [_simplify_ast(term) for term in pos]
        neg = [_simplify_ast(term) for term in neg]
        
        # We collect and sum the constant values
        values = [term.value for term in pos if isinstance(term, Constant)] +\
                [-term.value for term in neg if isinstance(term, Constant)] 
        value = sum(values)
        # print("value", value)
        # Remove the constants from our pos and neg lists
        pos = [term for term in pos if not isinstance(term, Constant)]
        neg = [term for term in neg if not isinstance(term, Constant)]
        # print("pos", dump(pos[0]))
        # print("neg", dump(neg[0]))
        new_pos, new_neg = [], []
        for term in pos:
            if isinstance(term, UnaryOp):
                if isinstance(term, USub):
                    new_neg.append(term.operand)
            else:
                # print("entered here", term)
                new_pos.append(term)
        for term in neg:
            if isinstance(term, UnaryOp):
                if isinstance(term, USub):
                    new_pos.append(term.operand)
            else:
                # print("entered here2", term)
                new_neg.append(term)
        pos, neg = new_pos, new_neg
        # print("pos1", pos)
        # print("neg1", neg)
        # Append the constant value sum to pos or neg
        if value > 0:
            pos.append(Constant(value=value))
        elif value < 0:
            neg.append(Constant(value=abs(value)))
        # print("pos3", pos)
        # print("neg3", neg)
        # Count the number of occurances of each term.
        term_counts = [(term, get_count_from_ast(pos, term) - get_count_from_ast(neg, term)) for term in
                       pos + neg]
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        # print(term_counts)
        term_counts = dict([(AST.ast2str(term), (term, count)) for term, count in 
                            term_counts])
        # print("after",term_counts)
        # We find the first term with non-zero count.
        ii = 0
        for ii, term in enumerate(pos+neg):
            ast_out, count = term_counts[AST.ast2str(term)]
            if count != 0:
                break
        else:
            # We get here if we don't break out of the loop, implying that
            #  all our terms had count of 0
            return _ZERO
        # print("first", dump(ast_out))
        term_counts[AST.ast2str(term)] = (ast_out, 0)
        if abs(count) != 1:
            ast_out = BinOp(left=Constant(value=abs(count)), op=Mult(), right=ast_out)
        if count < 0:
            ast_out = UnaryOp(op=USub(), operand=ast_out)

        # And add in all the rest
        for term in (pos+neg)[ii:]:
            term, count = term_counts[AST.ast2str(term)]
            term_counts[AST.ast2str(term)] = (term, 0)
            if abs(count) != 1:
                term = BinOp(left=Constant(value=abs(count)), op=Mult(), right=term)
            if count > 0:
                ast_out = BinOp(left=ast_out, op=Add(), right=term)
            elif count < 0:
                ast_out = BinOp(left=ast_out, op=Sub(), right=term)
        # print("ast_out", dump(ast_out))
        return ast_out
    
    elif isinstance(ast, BinOp) and (isinstance(ast.op, Mult) or isinstance(ast.op, Div)):
        # print("entered########################################", unparse(ast))
        # We collect numerator and denominator terms and simplify each of them
        num, denom = [], []
        AST._collect_num_denom(ast, num, denom)
        # print("hhhhhhhhhhhhhhhhhhhhhhhhhhhhh")
        # print("num", num)
        # print("denom", denom)
        # print("first", dump(num[0]), dump(num[1]))
        num = [_simplify_ast(term) for term in num]
        denom = [_simplify_ast(term) for term in denom]
        # print("num1", num)
        # print("denom2", denom)
        # We collect and sum the constant values
        values = [term.value for term in num if isinstance(term, Constant)] +\
                [1./term.value for term in denom if isinstance(term, Constant)] 
        # This takes the product of all our values
        # print("values before functools", values)
        value = functools.reduce(operator.mul, values + [1])
        # print("value", value)
        # If our value is 0, the expression is 0
        if not value:
            return _ZERO
        # print("num2", num)
        # print("denom2", denom)
        # Remove the constants from our pos and neg lists
        num = [term for term in num if not isinstance(term, Constant)]
        denom = [term for term in denom if not isinstance(term, Constant)]
        # print("after", dump(num[0]), dump(num[1]))
        # print("denom", denom)
        # Here we count all the negative (UnarySub) elements of our expression.
        # We also remove the UnarySubs from their arguments. We'll correct
        #  for it at the end.
        num_neg = 0
        for list_of_terms in [num, denom]:
            for ii, term in enumerate(list_of_terms):
                if isinstance(term, USub):
                    list_of_terms[ii] = term.operand
                    num_neg += 1

        # Append the constant value sum to pos or neg
        if abs(value) != 1:
            num.append(Constant(value=abs(value)))
        if value < 0:
            num_neg += 1

        make_neg = num_neg % 2
        # for i in num:
        #     print("testttttttttttttttttttttttttttttttttttttttttttttttttttt", dump(i))
        # Count the number of occurances of each term.
        term_counts = [(term, get_count_from_ast(num, term) - get_count_from_ast(denom, term)) for term in
                       num + denom]
        # print("qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq")
        # print(term_counts)
        # print("term counts", term_counts)
        # Tricky: We use the str(term) as the key for the dictionary to ensure
        #         that each entry represents a unique term. We also drop terms
        #         that have a total count of 0.
        term_counts = dict([(AST.ast2str(term), (term, count)) for term, count in 
                            term_counts])

        nums, denoms = [], []
        # We walk through terms in num+denom in order, so we rearrange a little
        #  as possible.
        for term in num+denom:
            term, count = term_counts[AST.ast2str(term)]
            # Once a term has been done, we set its term_counts to 0, so it
            #  doesn't get done again.
            term_counts[AST.ast2str(term)] = (term, 0)
            if abs(count) > 1:
                term = BinOp(left=term, op=Pow(), right=Constant(value=abs(count)))
                # term = Pow((term, Constant(abs(count))))
            if count > 0:
                nums.append(term)
            elif count < 0:
                denoms.append(term)

        # We return the product of the numerator terms over the product of the
        #  denominator terms
        out = AST._make_product(nums)
        if denoms:
            denom = AST._make_product(denoms)
            # out = Div((out, denom))
            out = BinOp(left=out, op=Div(), right=denom)

        if make_neg:
            out = UnaryOp(op=USub(), operand=out)
        # print("here-----------***********************", dump(out))
        return out
    elif isinstance(ast, BinOp) and isinstance(ast.op, Pow):
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
        elif isinstance(base, Constant) and\
                isinstance(power, Constant):
            return Constant(value=base.value**power.value)
        # Getting here implies that no simplifications are possible, so just
        #  return with simplified arguments
        # print("in power &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
        # print(dump(BinOp(left=base, op=Pow(), right=power)))
        return BinOp(left=base, op=Pow(), right=power)
    
    elif isinstance(ast, UnaryOp) and isinstance(ast.op, USub):
        simple_expr = _simplify_ast(ast.operand)
        if isinstance(simple_expr, UnaryOp) and isinstance(simple_expr.op, USub):
            # Case --x
            return _simplify_ast(simple_expr.operand)
        elif isinstance(simple_expr, Constant):
            if simple_expr.value == 0:
                return Constant(value=0)
            else:
                return Constant(value=-simple_expr.value)
        else:
            return UnaryOp(op=USub(), operand=simple_expr)
    elif isinstance(ast, UnaryOp) and isinstance(ast.op, UAdd):
        simple_expr = _simplify_ast(ast.operand)
        return simple_expr
    elif isinstance(ast, list):
        simple_list = [_simplify_ast(elem) for elem in ast]
        return simple_list
    elif isinstance(ast, tuple):
        return tuple(_simplify_ast(list(ast)))
    elif ast.__class__ in AST._node_attrs:
        # print("}}}}}}}}}}}}}}}}}}}}}}}}}}}}")
        # Handle node types with no special cases.
        for attr_name in AST._node_attrs[ast.__class__]:
            # print("hhhhheeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee")
            # print(attr_name)
            attr = getattr(ast, attr_name)
            # print("attribute name", attr)
            if isinstance(attr, list):
                for ii, elem in enumerate(attr):
                    # print("before -----------------------------------------", ii, elem)
                    attr[ii] = _simplify_ast(elem)
                    # print("after--------------------------------------", attr[ii])
            else:
                setattr(ast, attr_name, _simplify_ast(attr))
        return ast
    else:
        return ast
