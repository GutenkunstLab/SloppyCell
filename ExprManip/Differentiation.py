from compiler.ast import *
import copy
import cPickle
import logging
logger = logging.getLogger('ExprManip.Differentiation')

import AST
from AST import strip_parse
import Simplify
import Substitution

_ZERO = Const(0)
_ONE = Const(1)

__deriv_saved = {}
def load_derivs(filename):
    global __deriv_saved
    try:
        __deriv_saved = cPickle.load(file(filename, 'rb'))
    except IOError:
        pass

def save_derivs(filename):
    f = file(filename, 'wb')
    cPickle.dump(__deriv_saved, f, 2)
    f.close()

def diff_expr(expr, wrt):
    """
    Return the derivative of the expression with respect to a given variable.
    """
    logger.debug('Taking derivative of %s wrt %s' % (expr, wrt))
    key = '%s__derivWRT__%s' % (expr, wrt)
    if __deriv_saved.has_key(key):
        deriv = __deriv_saved[key]
        logger.debug('Found saved result %s.' % deriv)
        return deriv

    ast = AST.strip_parse(expr)
    deriv = _diff_ast(ast, wrt)
    deriv = Simplify._simplify_ast(deriv)
    deriv  = AST.ast2str(deriv)
    __deriv_saved[key] = deriv
    logger.debug('Computed result %s.' % deriv)
    return deriv

# This dictionary stores how to differentiate various functions. The keys are 
#  (function name, # of arguments). The values are tuples of strings which 
#  give the partial derivive wrt each argument. It is important to use arg# to
#  denote the arguments.
_KNOWN_FUNCS = {('acos', 1): ('1/sqrt(1-arg0**2)',),
                ('asin', 1): ('1/sqrt(1-arg0**2)',),
                ('atan', 1): ('1/(1+arg0**2)',),
                ('cos', 1): ('-sin(arg0)',),
                ('cosh', 1): ('sinh(arg0)',),
                ('exp', 1): ('exp(arg0)',),
                ('log', 1): ('1/arg0',),
                ('log10', 1): ('1/(log(10)*arg0)',),
                ('sin', 1): ('cos(arg0)',),
                ('sinh', 1): ('cosh(arg0)',),
                ('sqrt', 1): ('1/(2*sqrt(arg0))',),
                ('tan', 1): ('1/cos(arg0)**2',),
                ('tanh', 1): ('1/cosh(arg0)**2',),
                ('pow', 2): ('arg1 * arg0**(arg1-1)', 'log(arg0) * arg0**arg1')
                }
for key, terms in _KNOWN_FUNCS.items():
    _KNOWN_FUNCS[key] = [strip_parse(term) for term in terms]

def _diff_ast(ast, wrt):
    """
    Return an AST that is the derivative of ast with respect the variable with
    name 'wrt'.
    """

    # For now, the strategy is to return the most general forms, and let 
    #  the simplifier take care of the special cases.
    if isinstance(ast, Name):
        if ast.name == wrt:
            return _ONE
        else:
            return _ZERO
    elif isinstance(ast, Const):
        return _ZERO
    elif isinstance(ast, Add) or isinstance(ast, Sub):
        # Just take the derivative of the arguments. The call to ast.__class__
        #  lets us use the same code from Add and Sub.
        return ast.__class__((_diff_ast(ast.left, wrt), 
                              _diff_ast(ast.right, wrt)))
    elif isinstance(ast, Mul) or isinstance(ast, Div):
        # Collect all the numerators and denominators together
        nums, denoms = [], []
        AST._collect_num_denom(ast, nums, denoms)

        # Collect the numerator terms into a single AST
        num = AST._make_product(nums)
        # Take the derivative of the numerator terms as a product
        num_d = _product_deriv(nums, wrt)
        if not denoms:
            # If there is no denominator
            return num_d

        denom = AST._make_product(denoms)
        denom_d = _product_deriv(denoms, wrt)

        # Derivative of x/y is x'/y + -x*y'/y**2
        term1 = Div((num_d, denom))
        term2 = Div((Mul((UnarySub(num), denom_d)), Power((denom, Const(2)))))
        return Add((term1, term2))

    elif isinstance(ast, Power):
        # Use the derivative of the 'pow' function
        ast =  CallFunc(Name('pow'), [ast.left, ast.right])
        return _diff_ast(ast, wrt)

    elif isinstance(ast, CallFunc):
        func_name = AST.ast2str(ast.node)
        args = ast.args
        args_d = [_diff_ast(arg, wrt) for arg in args]

        if _KNOWN_FUNCS.has_key((func_name, len(args))):
            form = copy.deepcopy(_KNOWN_FUNCS[(func_name, len(args))])
        else:
            # If this isn't a known function, our form is
            #  (f_0(args), f_1(args), ...)
            args_expr = [Name('arg%i' % ii) for ii in range(len(args))]
            form = [CallFunc(Name('%s_%i' % (func_name, ii)), args_expr) for
                    ii in range(len(args))]

        # We build up the terms in our derivative
        #  f_0(x,y)*x' + f_1(x,y)*y', etc.
        outs = []
        for arg_d, arg_form_d in zip(args_d, form):
            # We skip arguments with 0 derivative
            if arg_d == _ZERO:
                continue
            for ii, arg in enumerate(args):
                Substitution._sub_subtrees_for_vars(arg_form_d, 
                                                    {'arg%i'%ii:arg})
            outs.append(Mul((arg_form_d, arg_d)))

        # If all arguments had zero deriviative
        if not outs:
            return _ZERO
        else:
            # We add up all our terms
            ret = outs[0]
            for term in outs[1:]:
                ret = Add((ret, term))
            return ret

    elif isinstance(ast, UnarySub):
        return UnarySub(_diff_ast(ast.expr, wrt))

    elif isinstance(ast, UnaryAdd):
        return UnaryAdd(_diff_ast(ast.expr, wrt))

def _product_deriv(terms, wrt):
    """
    Return an AST expressing the derivative of the product of all the terms.
    """
    if len(terms) == 1:
        return _diff_ast(terms[0], wrt)
    deriv_terms = []
    for ii, term in enumerate(terms):
        term_d = _diff_ast(term, wrt)
        other_terms = terms[:ii] + terms[ii+1:]
        deriv_terms.append(AST._make_product(other_terms + 
                                                          [term_d]))
    sum = deriv_terms[0]
    for term in deriv_terms[1:]:
        sum = Add((term, sum))

    return sum
