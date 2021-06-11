from ast import *

TINY = 1e-12
        
# def _node_equal(self, other):
#     """
#     Return whether self and other represent the same expressions.

#     Unfortunately, the Node class in Python 2.3 doesn't define ==, so
#     we need to write our own.
#     """
#     # We're not equal if other isn't a Node, or if other is a different class.
#     if not isinstance(other, Node) or not isinstance(other, self.__class__):
#         return False
#     # Loop through all children, checking whether they are equal
#     for self_child, other_child in zip(self.getChildren(), other.getChildren()):
#         if not self_child == other_child:
#             return False
#     # If we get here, our two nodes much be equal
#     return True
# # Add our new equality tester to the Node class.
# Node = AST
# Node.__eq1__ = _node_equal
# class Visitor(NodeVisitor):
#     def visit(self, node):
#        if isinstance(node, self.whitelist):
#             return super().visit(node)

#     whitelist = (Name,Constant, Call, Subscript, Slice, Slice,Pow,USub,
#                  UAdd,Mult,Div,Sub,Add,Compare, Not, And, Or, Expr)
    
def strip_parse(expr):
    """
    Return an abstract syntax tree (AST) for an expression.

    This removes the enclosing cruft from a call to compiler.parse(expr)
    """
    # The .strip() ignores leading and trailing whitespace, which would 
    #  otherwise give syntax errors.
    tree = parse(str(expr).strip())
    return  tree

# This defines the order of operations for the various node types, to determine
#  whether or not parentheses are necessary.
_OP_ORDER = {Name: 0,
             Call: 0,
             Subscript: 0,
             Slice: 0,
             Slice: 0,
             Pow: 3,
             USub: 4,
             UAdd: 4,
             Mult: 5,
             Div: 5,
             Sub: 10,
             Add: 10,
             Compare: 11,
             Not: 11,
             And: 11,
             Or: 11,
             Expr: 100
            }

# This is just an instance of Discard to use for the default
_FARTHEST_OUT = Expr(None)

# These are the attributes of each node type that are other nodes.
_node_attrs = {Name: (),
               Add: ('left', 'right'),
               Sub: ('left', 'right'),
               Mult: ('left', 'right'),
               Div: ('left', 'right'),
               Call: ('args',),
               Pow: ('left', 'right'),
               USub: ('expr',),
               UAdd: ('expr',),
               Slice: ('lower', 'upper'),
               Slice: ('nodes',),
               Subscript: ('subs',),
               Compare: ('expr', 'ops'),
               Not: ('expr',),
               Or: ('nodes',),
               And: ('nodes',),
               }



# class Visitor(NodeTransformer):
#     def visit_BinOp(self, node):
#         node.left = self.visit(node.left)
#         node.right = self.visit(node.right)

#         if isinstance(node.op, Add):
#             node = '%s + %s' % (node.left,
#                            node.right)
#         return node

#     def visit_Sub(self, node):
#         print(node)
#         self.generic_visit(node)
        
#     def visit_Constant(self, node):
#         out = str(node.value)
#         self.generic_visit(node)

    


    

    

    

    
   

   

    
    



SOURCE = """
def hello(msg):
    a = 21 * 2
    print(msg, a)
"""



def ast2str(node, outer = _FARTHEST_OUT , adjust = 0):
    """
    Return the string representation of an AST.

    outer: The AST's 'parent' node, used to determine whether or not to 
        enclose the result in parentheses. The default of _FARTHEST_OUT will
        never enclose the result in parentheses.

    adjust: A numerical value to adjust the priority of this ast for
        particular cases. For example, the denominator of a '/' needs 
        parentheses in more cases than does the numerator.
    """
    return unparse(node)
    if isinstance(node, Name):
        out = node.name
    elif isinstance(node, Constant):
        out = str(node.value)
    elif isinstance(node, Add):
        out = '%s + %s' % (ast2str(node.left, node),
                           ast2str(node.right, node))
    elif isinstance(node, Sub):
        out = '%s - %s' % (ast2str(node.left, node),
                           ast2str(node.right, node, adjust = TINY))
    elif isinstance(node, Mult):
        out = '%s*%s' % (ast2str(node.left, node),
                           ast2str(node.right, node))
    elif isinstance(node, Div):
        # The adjust ensures proper parentheses for x/(y*z)
        out = '%s/%s' % (ast2str(node.left, node),
                           ast2str(node.right, node, adjust = TINY))
    elif isinstance(node, Pow):
        # The adjust ensures proper parentheses for (x**y)**z
        out = '%s**%s' % (ast2str(node.left, node, adjust = TINY),
                          ast2str(node.right, node))
    elif isinstance(node, USub):
        out = '-%s' % ast2str(node.expr, node)
    elif isinstance(node, UAdd):
        out = '+%s' % ast2str(node.expr, node)
    elif isinstance(node, Call):
        args = [ast2str(arg) for arg in node.args]
        out = '%s(%s)' % (ast2str(node.node), ', '.join(args))
    elif isinstance(node, Subscript):
        subs = [ast2str(sub) for sub in node.subs]
        out = '%s[%s]' % (ast2str(node.expr), ', '.join(subs))
    elif isinstance(node, Slice):
        out = '%s[%s:%s]' % (ast2str(node.expr), ast2str(node.lower), 
                             ast2str(node.upper))
    elif isinstance(node, Slice):
        nodes = [ast2str(node) for node in node.nodes]
        out = ':'.join(nodes)
    elif isinstance(node, Compare):
        expr = ast2str(node.expr, node, adjust=6+TINY)
        out_l = [expr]
        for op, val in node.ops:
            out_l.append(op)
            out_l.append(ast2str(val, node, adjust=6+TINY))
        out = ' '.join(out_l)
    elif isinstance(node, And):
        nodes = [ast2str(node, node, adjust=TINY) for node in node.nodes]
        out = ' and '.join(nodes)
    elif isinstance(node, Or):
        nodes = [ast2str(node, node, adjust=TINY) for node in node.nodes]
        out = ' or '.join(nodes)
    elif isinstance(node, Not):
        out = 'not %s' % ast2str(node.expr, node, adjust=TINY)
    # Ensure parentheses by checking the _OP_ORDER of the outer and inner ASTs
    if _need_parens(outer, node, adjust):
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
    print(outer, inner, adjust)
    return _OP_ORDER[outer.__class__] >= _OP_ORDER[inner.__class__] + adjust

def _collect_num_denom(node, nums, denoms):
    """
    Append to nums and denoms, respectively, the nodes in the numerator and 
    denominator of an AST.
    """
    if not (isinstance(node, Mult) or isinstance(node, Div)):
        # If ast is not multiplication or division, just put it in nums.
        nums.append(node)
        return

    if isinstance(node.left, Div) or isinstance(node.left, Mult):
        # If the left argument is a multiplication or division, descend into
        #  it, otherwise it is in the numerator.
        _collect_num_denom(node.left, nums, denoms)
    else:
        nums.append(node.left)

    if isinstance(node.right, Div) or isinstance(node.right, Mult):
        # If the left argument is a multiplication or division, descend into
        #  it, otherwise it is in the denominator.
        if isinstance(node, Mult):
            _collect_num_denom(node.left, nums, denoms)
        elif isinstance(node, Div):
            # Note that when we descend into the denominator of a Div, we want 
            #  to swap our nums and denoms lists
            _collect_num_denom(node.right, denoms, nums)
    else:
        if isinstance(node, Mult):
            nums.append(node.right)
        elif isinstance(node, Div):
            denoms.append(node.right)

def _collect_pos_neg(node, poss, negs):
    """
    Append to poss and negs, respectively, the nodes in AST with positive and 
    negative factors from a addition/subtraction chain.
    """
    # 
    # This code is almost duplicated from _collect_num_denom. 
    # The additional twist is handling UnarySubs.
    #
    if not (isinstance(node, Add) or isinstance(node, Sub)):
        poss.append(node)
        return

    if isinstance(node.left, Sub) or isinstance(node.left, Add):
        _collect_pos_neg(node.left, poss, negs)
    else:
        poss.append(node.left)

    if isinstance(node.right, Sub) or isinstance(node.right, Add):
        if isinstance(node, Add):
            _collect_pos_neg(node.right, poss, negs)
        elif isinstance(node, Sub):
            _collect_pos_neg(node.right, negs, poss)
    else:
        if isinstance(node, Add):
            poss.append(node.right)
        elif isinstance(node, Sub):
            negs.append(node.right)

def _make_product(terms):
    """
    Return an AST expressing the product of all the terms.
    """
    if terms:
        product = terms[0]
        for term in terms[1:]:
            product = Mult((product, term))
        return product 
    else:
        return Constant(1)

def recurse_down_tree(ast, func, args=()):
    if isinstance(ast, list):
        for ii, elem in enumerate(ast):
            ast[ii] = func(elem, *args)
    elif isinstance(ast, tuple):
        ast = tuple(func(list(ast), *args))
    elif ast.__class__ in _node_attrs:
        for attr_name in _node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            attr_mod = func(attr, *args)
            setattr(ast, attr_name, attr_mod)

    return ast

