from compiler.ast import *
import sets

import AST

extract_vars_cache = {}

def extract_comps(expr):
    """
    Extract all comparisons from the expression.
    """
    comps_found = []
    _extract_comps_ast(AST.strip_parse(expr), comps_found)
    comps_found = [AST.ast2str(ast) for ast in comps_found]
    return sets.Set(comps_found)

def _extract_comps_ast(ast, comps_found):
    if isinstance(ast, Compare):
        comps_found.append(ast)
        _extract_comps_ast(ast.expr, comps_found)
        for op, elem in ast.ops:
            _extract_comps_ast(elem, comps_found)
    elif isinstance(ast, list) or isinstance(ast, tuple):
        for elem in ast:
            _extract_comps_ast(elem, comps_found)
    elif AST._node_attrs.has_key(ast.__class__):
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            _extract_comps_ast(attr, comps_found)

def extract_vars(expr):
    """
    Return a Set of the variables used in an expression.
    """
    try:
        return extract_vars_cache[expr]
    except KeyError:
        vars_found = []
        _extract_vars_ast(AST.strip_parse(expr), vars_found)
        vars_found = [AST.ast2str(ast) for ast in vars_found]
        result = sets.Set(vars_found)
        extract_vars_cache[expr] = result
        return result

def _extract_vars_ast(ast, vars_found):
    """
    Appends the asts of the variables used in ast to vars_found.
    """
    if isinstance(ast, Name):
        if ast.name not in ['True', 'False']:
            vars_found.append(ast)
    ast = AST.recurse_down_tree(ast, _extract_vars_ast, (vars_found,))
    return ast

def extract_funcs(expr):
    """
    Return a Set of the functions used in an expression.

    The elements of the Set are ('function_name', #arguments).
    """
    funcs_found = []
    _extract_funcs_ast(AST.strip_parse(expr), funcs_found)
    return sets.Set(funcs_found)

def _extract_funcs_ast(ast, funcs_found):
    """
    Append ('name', #arg) for each function used in the ast to funcs_found.
    """
    if isinstance(ast, CallFunc):
        funcs_found.append((AST.ast2str(ast.node), len(ast.args)))
        for node in ast.args:
            _extract_funcs_ast(node, funcs_found)
    ast = AST.recurse_down_tree(ast, _extract_funcs_ast, (funcs_found,))
    return ast
