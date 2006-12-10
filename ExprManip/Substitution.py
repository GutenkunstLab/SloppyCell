from compiler.ast import *
import copy

import AST
from AST import strip_parse, ast2str
import Simplify

def sub_for_comps(expr, mapping):
    """
    For each pair out_name:in_expr in mapping, the returned string has all
    occurences of the variable out_compe substituted by in_expr.
    """
    if len(mapping) == 0:
        return expr

    ast = strip_parse(expr)
    ast_mapping = {}
    for out_expr, in_expr in mapping.items():
        out_ast = strip_parse(out_expr)
        if not isinstance(out_ast, Compare):
            raise ValueError, 'Expression %s to substitute for is not a '\
                    'comparison.' % out_expr
        ast_mapping[ast2str(out_ast)] = strip_parse(in_expr)

    ast = _sub_subtrees_for_comps(ast, ast_mapping)
    return ast2str(ast)

def _sub_subtrees_for_comps(ast, ast_mappings):
    if isinstance(ast, Compare) and ast_mappings.has_key(ast2str(ast)):
        return ast_mappings[ast2str(ast)]
    elif isinstance(ast, list):
        ast = [_sub_subtrees_for_comps(elem, ast_mappings) for elem in ast]
    elif isinstance(ast, tuple):
        ast = _sub_subtrees_for_comps(list(ast), ast_mappings)
        ast = tuple(ast)
    elif AST._node_attrs.has_key(ast.__class__):
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            setattr(ast, attr_name, _sub_subtrees_for_comps(attr, ast_mappings))
    return ast

def sub_for_var(expr, out_name, in_expr):
    """
    Returns a string with all occurances of the variable out_name substituted by
    in_expr.
    
    Perhaps regular expressions could do this more simply...
    """
    return sub_for_vars(expr, {out_name:in_expr})

def sub_for_vars(expr, mapping):
    """
    For each pair out_name:in_expr in mapping, the returned string has all
    occurences of the variable out_name substituted by in_expr.
    """
    if len(mapping) == 0:
        return expr

    ast = strip_parse(expr)
    ast_mapping = {}
    for out_name, in_expr in mapping.items():
        out_ast = strip_parse(out_name)
        if not isinstance(out_ast, Name):
            raise ValueError, 'Expression %s to substitute for is not a '\
                    'variable name.' % out_name
        ast_mapping[str(out_ast.name)] = strip_parse(in_expr)

    ast = _sub_subtrees_for_vars(ast, ast_mapping)
    return ast2str(ast)

def _sub_subtrees_for_vars(ast, ast_mappings):
    """
    For each out_name, in_ast pair in mappings, substitute in_ast for all 
    occurances of the variable named out_name in ast
    """
    if isinstance(ast, Name) and ast_mappings.has_key(ast2str(ast)):
        return ast_mappings[ast2str(ast)]
    elif isinstance(ast, list):
        ast = [_sub_subtrees_for_vars(elem, ast_mappings) for elem in ast]
    elif isinstance(ast, tuple):
        ast = _sub_subtrees_for_vars(list(ast), ast_mappings)
        ast = tuple(ast)
    elif AST._node_attrs.has_key(ast.__class__):
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            setattr(ast, attr_name, _sub_subtrees_for_vars(attr, ast_mappings))
    return ast

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
        mapping = {}
        for var_name, arg_ast in zip(func_vars, ast.args):
            subbed_arg_ast = _sub_for_func_ast(arg_ast, func_name, func_vars, 
                                               func_expr_ast)
            mapping[var_name] = subbed_arg_ast
        _sub_subtrees_for_vars(working_ast, mapping)
        return working_ast
    # Else we walk through the attributes of our ast, checking whether they
    #  need to be substituted. We do this because we can't, in general,
    #  substitute the function in-place.
    elif isinstance(ast, list):
        for ii, elem in enumerate(ast):
            ast[ii] = _sub_for_func_ast(elem, func_name, func_vars, 
                                        func_expr_ast)
    elif isinstance(ast, tuple):
        ast = tuple(_sub_for_func_ast(list(ast), func_name, func_vars,
                                      func_expr_ast))
    elif AST._node_attrs.has_key(ast.__class__):
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            attr_mod = _sub_for_func_ast(attr, func_name, func_vars,
                                         func_expr_ast)
            setattr(ast, attr_name, attr_mod)
                    
    return ast
