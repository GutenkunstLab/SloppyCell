import compiler
from compiler.ast import *

TINY = 1e-12

def _node_equal(self, other):
    """
    Return whether self and other represent the same expressions.

    Unfortunately, the Node class in Python 2.3 doesn't define ==, so
    we need to write our own.
    """
    # We're not equal if other isn't a Node, or if other is a different class.
    if not isinstance(other, Node) or not isinstance(other, self.__class__):
        return False
    # Loop through all children, checking whether they are equal
    for self_child, other_child in zip(self.getChildren(), other.getChildren()):
        if not self_child == other_child:
            return False
    # If we get here, our two nodes much be equal
    return True
# Add our new equality tester to the Node class.
Node.__eq__ = _node_equal

def strip_parse(expr):
    """
    Return an abstract syntax tree (AST) for an expression.

    This removes the enclosing cruft from a call to compiler.parse(expr)
    """
    # The .strip() ignores leading and trailing whitespace, which would 
    #  otherwise give syntax errors.
    ast = compiler.parse(str(expr).strip())
    return ast.node.nodes[0].expr

# This defines the order of operations for the various node types, to determine
#  whether or not parentheses are necessary.
_OP_ORDER = {Name: 0,
             Const: 0,
             CallFunc: 0,
             Subscript: 0,
             Slice: 0,
             Sliceobj: 0,
             Power: 3,
             UnarySub: 4,
             UnaryAdd: 4,
             Mul: 5,
             Div: 5,
             Sub: 10,
             Add: 10,
             Compare: 11,
             Not: 11,
             And: 11,
             Or: 11,
             Discard: 100}

# This is just an instance of Discard to use for the default
_FARTHEST_OUT = Discard(None)

# These are the attributes of each node type that are other nodes.
_node_attrs = {Name: (),
               Const: (),
               Add: ('left', 'right'),
               Sub: ('left', 'right'),
               Mul: ('left', 'right'),
               Div: ('left', 'right'),
               CallFunc: ('args',),
               Power: ('left', 'right'),
               UnarySub: ('expr',),
               UnaryAdd: ('expr',),
               Slice: ('lower', 'upper'),
               Sliceobj: ('nodes',),
               Subscript: ('subs',),
               Compare: ('expr', 'ops'),
               Not: ('expr',),
               Or: ('nodes',),
               And: ('nodes',),
               }


def ast2str(ast, outer = _FARTHEST_OUT, adjust = 0):
    """
    Return the string representation of an AST.

    outer: The AST's 'parent' node, used to determine whether or not to 
        enclose the result in parentheses. The default of _FARTHEST_OUT will
        never enclose the result in parentheses.

    adjust: A numerical value to adjust the priority of this ast for
        particular cases. For example, the denominator of a '/' needs 
        parentheses in more cases than does the numerator.
    """
    if isinstance(ast, Name):
        out = ast.name
    elif isinstance(ast, Const):
        out = str(ast.value)
    elif isinstance(ast, Add):
        out = '%s + %s' % (ast2str(ast.left, ast),
                           ast2str(ast.right, ast))
    elif isinstance(ast, Sub):
        out = '%s - %s' % (ast2str(ast.left, ast),
                           ast2str(ast.right, ast, adjust = TINY))
    elif isinstance(ast, Mul):
        out = '%s*%s' % (ast2str(ast.left, ast),
                           ast2str(ast.right, ast))
    elif isinstance(ast, Div):
        # The adjust ensures proper parentheses for x/(y*z)
        out = '%s/%s' % (ast2str(ast.left, ast),
                           ast2str(ast.right, ast, adjust = TINY))
    elif isinstance(ast, Power):
        # The adjust ensures proper parentheses for (x**y)**z
        out = '%s**%s' % (ast2str(ast.left, ast, adjust = TINY),
                          ast2str(ast.right, ast))
    elif isinstance(ast, UnarySub):
        out = '-%s' % ast2str(ast.expr, ast)
    elif isinstance(ast, UnaryAdd):
        out = '+%s' % ast2str(ast.expr, ast)
    elif isinstance(ast, CallFunc):
        args = [ast2str(arg) for arg in ast.args]
        out = '%s(%s)' % (ast2str(ast.node), ', '.join(args))
    elif isinstance(ast, Subscript):
        subs = [ast2str(sub) for sub in ast.subs]
        out = '%s[%s]' % (ast2str(ast.expr), ', '.join(subs))
    elif isinstance(ast, Slice):
        out = '%s[%s:%s]' % (ast2str(ast.expr), ast2str(ast.lower), 
                             ast2str(ast.upper))
    elif isinstance(ast, Sliceobj):
        nodes = [ast2str(node) for node in ast.nodes]
        out = ':'.join(nodes)
    elif isinstance(ast, Compare):
        expr = ast2str(ast.expr, ast, adjust=6+TINY)
        out_l = [expr]
        for op, val in ast.ops:
            out_l.append(op)
            out_l.append(ast2str(val, ast, adjust=6+TINY))
        out = ' '.join(out_l)
    elif isinstance(ast, And):
        nodes = [ast2str(node, ast, adjust=TINY) for node in ast.nodes]
        out = ' and '.join(nodes)
    elif isinstance(ast, Or):
        nodes = [ast2str(node, ast, adjust=TINY) for node in ast.nodes]
        out = ' or '.join(nodes)
    elif isinstance(ast, Not):
        out = 'not %s' % ast2str(ast.expr, ast, adjust=TINY)

    # Ensure parentheses by checking the _OP_ORDER of the outer and inner ASTs
    if _need_parens(outer, ast, adjust):
        return out
    else:
        return '(%s)' % out

def _need_parens(outer, inner, adjust):
    """
    Return whether or not the inner AST needs parentheses when enclosed by the
    outer.

    adjust: A numerical value to adjust the priority of this ast for
        particular cases. For example, the denominator of a '/' needs 
        parentheses in more cases than does the numerator.
    """
    return _OP_ORDER[outer.__class__] >= _OP_ORDER[inner.__class__] + adjust

def _collect_num_denom(ast, nums, denoms):
    """
    Append to nums and denoms, respectively, the nodes in the numerator and 
    denominator of an AST.
    """
    if not (isinstance(ast, Mul) or isinstance(ast, Div)):
        # If ast is not multiplication or division, just put it in nums.
        nums.append(ast)
        return

    if isinstance(ast.left, Div) or isinstance(ast.left, Mul):
        # If the left argument is a multiplication or division, descend into
        #  it, otherwise it is in the numerator.
        _collect_num_denom(ast.left, nums, denoms)
    else:
        nums.append(ast.left)

    if isinstance(ast.right, Div) or isinstance(ast.right, Mul):
        # If the left argument is a multiplication or division, descend into
        #  it, otherwise it is in the denominator.
        if isinstance(ast, Mul):
            _collect_num_denom(ast.right, nums, denoms)
        elif isinstance(ast, Div):
            # Note that when we descend into the denominator of a Div, we want 
            #  to swap our nums and denoms lists
            _collect_num_denom(ast.right, denoms, nums)
    else:
        if isinstance(ast, Mul):
            nums.append(ast.right)
        elif isinstance(ast, Div):
            denoms.append(ast.right)

def _collect_pos_neg(ast, poss, negs):
    """
    Append to poss and negs, respectively, the nodes in AST with positive and 
    negative factors from a addition/subtraction chain.
    """
    # 
    # This code is almost duplicated from _collect_num_denom. 
    # The additional twist is handling UnarySubs.
    #
    if not (isinstance(ast, Add) or isinstance(ast, Sub)):
        poss.append(ast)
        return

    if isinstance(ast.left, Sub) or isinstance(ast.left, Add):
        _collect_pos_neg(ast.left, poss, negs)
    else:
        poss.append(ast.left)

    if isinstance(ast.right, Sub) or isinstance(ast.right, Add):
        if isinstance(ast, Add):
            _collect_pos_neg(ast.right, poss, negs)
        elif isinstance(ast, Sub):
            _collect_pos_neg(ast.right, negs, poss)
    else:
        if isinstance(ast, Add):
            poss.append(ast.right)
        elif isinstance(ast, Sub):
            negs.append(ast.right)

def _make_product(terms):
    """
    Return an AST expressing the product of all the terms.
    """
    if terms:
        product = terms[0]
        for term in terms[1:]:
            product = Mul((product, term))
        return product 
    else:
        return Const(1)

def recurse_down_tree(ast, func, args=()):
    if isinstance(ast, list):
        for ii, elem in enumerate(ast):
            ast[ii] = func(elem, *args)
    elif isinstance(ast, tuple):
        ast = tuple(func(list(ast), *args))
    elif _node_attrs.has_key(ast.__class__):
        for attr_name in _node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            attr_mod = func(attr, *args)
            setattr(ast, attr_name, attr_mod)

    return ast

