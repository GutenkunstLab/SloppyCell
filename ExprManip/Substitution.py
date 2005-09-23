from compiler.ast import *
import copy

import AST
from AST import strip_parse, ast2str
import Simplify

def sub_for_var(expr, out_name, in_expr):
    """
    Returns a string with all occurances of the variable out_name substituted by
    in_expr.
    
    Perhaps regular expressions could do this more simply...
    """
    ast = strip_parse(expr)
    out_ast = strip_parse(out_name)
    if out_ast.__class__ != Name:
        raise ValueError, 'Expression to substitute for is not a variable name.'
    out_name = out_ast.name
    in_ast = strip_parse(in_expr)

    if ast.__class__ == Name and ast2str(ast) == out_name:
        return ast2str(in_ast)
    else:
        _sub_subtree_for_var(ast, out_name, in_ast)
        return ast2str(ast)

def _sub_subtree_for_var(ast, out_name, in_ast):
    """
    Substitute in_ast for all occurances of the variable named out_name in ast
    """
    for attr_name in AST._node_attrs[ast.__class__]:
        attr = getattr(ast, attr_name)
        if isinstance(attr, list):
            for ii, elem in enumerate(attr):
                if elem.__class__ == Name and ast2str(elem) == out_name:
                    attr[ii] = in_ast
                else:
                    _sub_subtree_for_var(elem, out_name, in_ast)
        else:
            if attr.__class__ == Name and ast2str(attr) == out_name:
                setattr(ast, attr_name, in_ast)
            else:
                _sub_subtree_for_var(attr, out_name, in_ast)

def sub_for_func(expr, func_name, func_vars, func_expr):
    """
    Return a string with the function func_name substituted for it's exploded 
    form.
    
    func_name: The name of the function.
    func_vars: A sequence variables used by the function expression
    func_expr: The expression for the function.

    For example:
        If f(x, y, z) = sqrt(z)*x*y-z
        func_name = 'f'
        func_vars = ['x', 'y', 'z']
        func_expr = 'sqrt(z)*x*y-z'

    Maybe I don't really neeed this... I should check speed trade off of just
    using lambdas in our code.
    """
    ast = strip_parse(expr)
    func_name_ast = strip_parse(func_name)
    if not isinstance(func_name_ast, Name):
        raise ValueError, 'Function name is not a simple name.'
    func_name = func_name_ast.name

    func_vars_ast = [strip_parse(var) for var in func_vars]
    for var_ast in func_vars_ast:
        if not isinstance(var_ast, Name):
            raise ValueError, 'Function variable is not a simple name.'
    func_var_names = [getattr(var_ast, 'name') for var_ast in func_vars_ast]

    func_expr_ast = strip_parse(func_expr)
    ast = _sub_for_func_ast(ast, func_name, func_var_names, func_expr_ast)
    simple = Simplify._simplify_ast(ast)
    return ast2str(simple)

def _sub_for_func_ast(ast, func_name, func_vars, func_expr_ast):
    """
    Return and ast with the function func_name substituted out.
    """
    if isinstance(ast, CallFunc) and ast2str(ast.node) == func_name\
       and len(ast.args) == len(func_vars):
        # If our ast is the function we're looking for, we take the ast
        #  for the function expression, substitute for its arguments, and
        #  return
        working_ast = copy.deepcopy(func_expr_ast)
        for var_name, arg_ast in zip(func_vars, ast.args):
            subbed_arg_ast = _sub_for_func_ast(arg_ast, func_name, func_vars, 
                                               func_expr_ast)
            _sub_subtree_for_var(working_ast, var_name, subbed_arg_ast)
        return working_ast
    # Else we walk through the attributes of our ast, checking whether they
    #  need to be substituted. We do this because we can't, in general,
    #  substitute the function in-place.
    else:
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            if isinstance(attr, list):
                for ii, elem in enumerate(attr):
                    attr[ii] = _sub_for_func_ast(elem, func_name, func_vars, 
                                                 func_expr_ast)
            else:
                setattr(ast, attr_name, _sub_for_func_ast(attr, func_name, 
                                                      func_vars, func_expr_ast))
    return ast
