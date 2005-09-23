from compiler.ast import *
import sets

import AST

def extract_vars(expr):
    """
    Return a Set of the variables used in an expression.
    """
    vars_found = []
    _extract_vars_ast(AST.strip_parse(expr), vars_found)
    vars_found = [AST.ast2str(ast) for ast in vars_found]
    return sets.Set(vars_found)

def _extract_vars_ast(ast, vars_found):
    """
    Append the asts of the variables used in ast to vars_found.
    """
    if isinstance(ast, Name):
        vars_found.append(ast)
    for attr_name in AST._node_attrs[ast.__class__]:
        attr = getattr(ast, attr_name)
        if isinstance(attr, list):
            for elem in attr:
                _extract_vars_ast(elem, vars_found)
        else:
            _extract_vars_ast(attr, vars_found)

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
    else:
        for attr_name in AST._node_attrs[ast.__class__]:
            attr = getattr(ast, attr_name)
            if isinstance(attr, list):
                for elem in attr:
                    _extract_funcs_ast(elem, funcs_found)
            else:
                _extract_funcs_ast(attr, funcs_found)
