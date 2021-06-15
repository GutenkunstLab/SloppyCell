from ast import *

TINY = 1e-12
    
def strip_parse(expr):
    """
    Return an abstract syntax tree (AST) for an expression.

    This removes the enclosing cruft from a call to compiler.parse(expr)
    """
    # The .strip() ignores leading and trailing whitespace, which would 
    #  otherwise give syntax errors.
    tree = parse(str(expr).strip())
    return  tree.body[0].value


def ast2str(ast):
    """
    Return the string representation of an AST.
    """
    return unparse(ast)
    

def _need_parens(outer, inner, adjust):
    """
    Return whether or not the inner AST needs parentheses when enclosed by the
    outer.

    adjust: A numerical value to adjust the priority of this ast for
        particular cases. For example, the denominator of a '/' needs 
        parentheses in more cases than does the numerator.
    """
    print(outer, inner, adjust)
    return _OP_ORDER[outer.__class__] >= _OP_ORDER[inner.__class__] + adjust

def _collect_num_denom(ast, nums, denoms):
    """
    Append to nums and denoms, respectively, the nodes in the numerator and 
    denominator of an AST.
    """
    if not isinstance(ast, BinOp):
        nums.append(ast.value)
        return
    if (isinstance(ast.op, Div) or isinstance(ast.op, Mult)):
        if isinstance(ast.left, BinOp):
            _collect_num_denom(ast.left, nums, denoms)
        else:
            nums.append(ast.left.value)
            
        if isinstance(ast.right, BinOp):
            if isinstance(ast.op, Div):
                _collect_num_denom(ast.right, denoms, nums)
            elif isinstance(ast.op, Mult):
                _collect_num_denom(ast.left, nums, denoms)
        else:
            if isinstance(ast.op, Div):
                denoms.append(ast.right.value)
            elif isinstance(ast.op, Mult):
                nums.append(ast.left.value)

def _collect_pos_neg(ast, poss, negs):
    """
    Append to poss and negs, respectively, the nodes in AST with positive and 
    negative factors from a addition/subtraction chain.
    """
    # 
    # This code is almost duplicated from _collect_num_denom. 
    # The additional twist is handling UnarySubs.
    #
    
    if not isinstance(ast, BinOp):
        poss.append(ast.value)
        return
    if (isinstance(ast.op, Sub) or isinstance(ast.op, Add)):
        if isinstance(ast.left, BinOp):
            _collect_pos_neg(ast.left, poss, negs)
        else:
            poss.append(ast.left.value)
            
        if isinstance(ast.right, BinOp):
            if isinstance(ast.op, Add):
                _collect_pos_neg(ast.right, poss, negs)
            elif isinstance(ast.op, Sub):
                _collect_pos_neg(ast.right, negs, poss)
                
        else:
            if isinstance(ast.op, Add):
                poss.append(ast.right.value)
            elif isinstance(ast.op, Sub):
                negs.append(ast.right.value)
