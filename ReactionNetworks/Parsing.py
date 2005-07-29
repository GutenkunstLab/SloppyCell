import copy, sets

import symbolic

descenders = ['power', 'atom', 'arith_expr', 'term', 'arglist', 'trailer', 'factor']

def extractVariablesFromString(input):
    ast = symbolic.string2ast(input)
    return extractVariablesFromAST(ast)

def extractVariablesFromAST(ast, variables = None):
    if variables is None:
        variables = sets.Set()

    # If we're at a leaf, just add the name and return
    if ast[0] is 'NAME':
        variables.add(ast[1])
        return variables

    # If we're in x**y
    if ast[0] is 'power' and ast[-2] == ['DOUBLESTAR', '**']:
        # Extract from x 
        extractVariablesFromAST(ast[:-2], variables)
        # Extract from y
        extractVariablesFromAST(ast[-1], variables)
        return variables

    # If we're in f(x)
    if ast[0] is 'power' and ast[-1][0] is 'trailer'\
       and ast[-1][1] == ['LPAR', '(']:
        # Extract from x
        extractVariablesFromAST(ast[-1][2], variables)
        return variables

    # If we're in x.y.z
    if ast[0] is 'power' and ast[-1][0] is 'trailer'\
       and ast[-1][1] == ['DOT', '.']:
        variables.add(symbolic.ast2string(ast))
        return variables

    # Else run through all our terms
    for term in ast:
        # Short-cut the recursion here
        if term[0] == 'NAME':
            variables.add(term[1])
        elif term[0] in descenders:
            extractVariablesFromAST(term, variables)

    return variables

def extractFunctionsFromString(input):
    ast = symbolic.string2ast(input)
    return extractFunctionsFromAST(ast)

def extractFunctionsFromAST(ast, functions = None):
    if functions is None:
        functions = sets.Set()

    # If we're in x**y
    if ast[0] is 'power' and ast[-2] == ['DOUBLESTAR', '**']:
        # Extract from x 
        extractFunctionsFromAST(ast[:-2], functions)
        # Extract from y
        extractFunctionsFromAST(ast[-1], functions)
        return functions

    # If we're in f(x)
    if ast[0] is 'power' and ast[-1][0] is 'trailer' and\
       ast[-1][1] == ['LPAR', '(']:
        # Add f
        # Extract the arguments to the function, discarding commas
        if ast[-1][2][0] == 'arglist':
            argsASTlist = ast[-1][2][1::2]
        else:
            argsASTlist = [ast[-1][2]]
        functions.add((symbolic.ast2string(ast[:-1]), len(argsASTlist)))

        # Extract from x
        extractFunctionsFromAST(ast[-1][2], functions)
        return functions

    for term in ast:
        if term[0] in descenders:
            extractFunctionsFromAST(term, functions)

    return functions

def substituteVariableNamesInString(input, old, new):
    ast = symbolic.string2ast(input)
    oldAST = symbolic.string2ast(old)
    newAST = symbolic.string2ast(new)
    substituteSubTreeInAST(ast, oldAST, newAST)
    return symbolic.ast2string(ast)

def substituteSubTreeInAST(ast, old, new):
    if ast == old:
        ast[:] = new
        return

    for ii, term in enumerate(ast):
        if term == old:
            ast[ii] = new
        elif term[0] in descenders:
            substituteSubTreeInAST(term, old, new)

def substituteFunctionIntoString(input, id, variables, expression):
    ast = symbolic.string2ast(input)
    idAST = symbolic.string2ast(id)
    if idAST[0] != 'power':
        idAST = ['power', idAST]

    variablesAST = map(symbolic.string2ast, variables)
    expressionAST = symbolic.string2ast(expression)
    ast = substituteFunctionIntoAST(ast, idAST, variablesAST, expressionAST)

    return symbolic.ast2string(ast)

def substituteFunctionIntoAST(ast, idAST, variablesAST, expressionAST):
    # If we're in x**y
    if ast[0] is 'power' and ast[-2] == ['DOUBLESTAR', '**']:
        # Subsitute into x 
        if ast[0] in descenders:
            ast[:-2] = substituteFunctionIntoAST(ast[:-2], idAST, variablesAST, 
                                                 expressionAST)
        # Substitute into y
        if ast[-1][0] in descenders:
            ast[-1] = substituteFunctionIntoAST(ast[-1], idAST, variablesAST, 
                                                expressionAST)

    # If we're in f(x) and we have the proper name
    if ast[0] is 'power' and ast[-1][0] is 'trailer' and\
       ast[-1][1] == ['LPAR', '('] and ast[:-1] == idAST:
        # Extract the arguments to the function, discarding commas
        if ast[-1][2][0] == 'arglist':
            argsASTlist = ast[-1][2][1::2]
        else:
            argsASTlist = [ast[-1][2]]

        if len(argsASTlist) == len(variablesAST):
            newExpressionAST = copy.deepcopy(expressionAST)
            for old, new in zip(variablesAST, argsASTlist):
                # If our argument is more than just a single number or variable 
                #  wrap it in parentheses
                if new[0] not in ['NAME', 'NUMBER']:
                    new = ['atom', ['LPAR', '('], new, ['RPAR', ')']]

                substituteSubTreeInAST(newExpressionAST, old, new)
                # Descend into the argument to see whether we need to substitute
                # there
                if new in descenders:
                    substituteFunctionIntoAST(new, idAST, variablesAST, 
                                              expressionAST)

            # Wrap our function up in parentheses to preserve order of 
            #  operations
            return ['atom', ['LPAR', '('], newExpressionAST, ['RPAR', ')']]

    for ii, term in enumerate(ast):
        if term[0] in descenders:
            ast[ii] = substituteFunctionIntoAST(term, idAST, variablesAST,
                                                expressionAST)

    return ast
#
# Tests
#

assert extractVariablesFromString('x') == sets.Set(['x'])
assert extractVariablesFromString('x + y') == sets.Set(['x', 'y'])
assert extractVariablesFromString('x * y') == sets.Set(['x', 'y'])
assert extractVariablesFromString('x / y') == sets.Set(['x', 'y'])
assert extractVariablesFromString('x**2') == sets.Set(['x'])
assert extractVariablesFromString('x**y') == sets.Set(['x', 'y'])
assert extractVariablesFromString('f(x)') == sets.Set(['x'])
assert extractVariablesFromString('f(x, y)') == sets.Set(['x', 'y'])
assert extractVariablesFromString('f(x + z/x.y, y)')\
        == sets.Set(['x', 'y', 'z', 'x.y'])
assert extractVariablesFromString('x**2 + f(y)') == sets.Set(['x', 'y'])
assert extractVariablesFromString('x + y**2 + z**(a + b) + z**f.g.h(a + b + sqrt(c))') == sets.Set(['x', 'y', 'z', 'a', 'b', 'c'])
assert extractVariablesFromString('g(a)**2') == sets.Set(['a'])

assert extractFunctionsFromString('f.g(x)') == sets.Set([('f.g', 1)])
assert extractFunctionsFromString('f(g(x))') == sets.Set([('f', 1), ('g', 1)])
assert extractFunctionsFromString('f(g(x)/2, a, b(x, y))')\
        == sets.Set([('f', 3), ('g', 1), ('b', 2)])
assert extractFunctionsFromString('a + g(x)') == sets.Set([('g', 1)])
assert extractFunctionsFromString('a + g(x) + (a + f(x))')\
        == sets.Set([('f', 1), ('g', 1)])
assert extractFunctionsFromString('f(g(1 + h.i(a))/2) + j(q**k(x)) + q.r[1] + s[1]') == sets.Set([('f', 1), ('g', 1), ('h.i', 1), ('j', 1), ('k', 1)])
assert extractFunctionsFromString('(g(x) + 2)**f(y)')\
        == sets.Set([('f', 1), ('g', 1)])
assert extractFunctionsFromString('g(x)**f(y)')\
        == sets.Set([('f', 1), ('g', 1)])

assert substituteVariableNamesInString('f(b, a) + a(q) + a.b(x) + a**x + y**a', 'a', 'c') == 'f(b,c)+c(q)+c.b(x)+c**x+y**c'
assert substituteVariableNamesInString('f.b(x + a.b)', 'a.b', 'c') == 'f.b(x+c)'
assert substituteVariableNamesInString('x', 'x', 'c') == 'c'

assert substituteFunctionIntoString('g(a + b, c)', 'g', ['x', 'y'], 'x + y**2')\
        == '((a+b)+c**2)'
assert substituteFunctionIntoString('f(g(a + b))', 'f', ['x'], 'h(x + 2)')\
        == '(h((g(a+b))+2))'
assert substituteFunctionIntoString('f.g.h(c)', 'f.g', ['x'], 'x + 2')\
        == 'f.g.h(c)'
assert substituteFunctionIntoString('f.g(c)', 'f.g', ['x'], 'x + 2') == '(c+2)'
assert substituteFunctionIntoString('f.g(c)**2', 'f.g', ['x'], 'x + 2')\
        == '(c+2)**2'
