from __future__ import nested_scopes
# Levenberg Marquardt minimization routines
"""
fmin_lm  : standard Levenberg Marquardt
fmin_lmNoJ : Levenberg Marquardt using a cost function instead of 
             a residual function and a gradient/J^tJ pair instead
             of the derivative of the residual function. Useful
             in problems where the number of residuals is very large.
fmin_lm_scale : scale invariant Levenberg Marquardt

"""
import scipy
from scipy import absolute, sqrt, asarray, zeros, mat, transpose, ones, dot, sum
import scipy.linalg
import copy
import SloppyCell.Utility
save = SloppyCell.Utility.save # module that provides pickled save

import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

abs = absolute
_epsilon = sqrt(scipy.finfo(scipy.float_).eps)

        
def approx_fprime(xk,f,epsilon,*args):
    f0 = apply(f,(xk,)+args)
    grad = scipy.zeros((len(xk),),scipy.float_)
    ei = scipy.zeros((len(xk),),scipy.float_)
    for k in range(len(xk)):
        ei[k] = epsilon
        grad[k] = (apply(f,(xk+ei,)+args) - f0)/epsilon
        ei[k] = 0.0
    return grad

def approx_fprime1(xk,f,epsilon,*args):
    """ centred difference formula to approximate fprime """
    #f0 = apply(f,(xk,)+args)
    grad = scipy.zeros((len(xk),),scipy.float_)
    ei = scipy.zeros((len(xk),),scipy.float_)
    epsilon = (epsilon**2.0)**(1.0/3.0)  # should be macheps^(1/3)
    for k in range(len(xk)):
        ei[k] = epsilon
        grad[k] = (apply(f,(xk+ei,)+args) - apply(f,(xk-ei,)+args))/(2.0*epsilon)
        ei[k] = 0.0
    return grad

def approx_fprime2(xk,f,epsilon,*args):
    """ centred difference formula to approximate the jacobian, given the residual 
    function """
    #f0 = apply(f,(xk,)+args)
    grad = scipy.zeros((len(xk),),scipy.float_)
    ei = scipy.zeros((len(xk),),scipy.float_)
    epsilon = (epsilon**2.0)**(1.0/3.0)  # should be macheps^(1/3)
    ei[0] = epsilon
    resminus = asarray(apply(f,(xk-ei,)+args))
    resplus = asarray(apply(f,(xk+ei,)+args))
    m = len(resminus)
    jac = scipy.zeros((m,len(xk)),scipy.float_)
    jac[:,0] = (resplus-resminus)/(2.0*epsilon)
    ei[0] = 0.0
    for k in range(1,len(xk)):
        ei[k] = epsilon
        resplus = asarray(apply(f,(xk+ei,)+args))
        resminus = asarray(apply(f,(xk-ei,)+args))
        jac[:,k] = (resplus-resminus)/(2.0*epsilon)
    #jac[k,:] = mat(transpose(mat(apply(f,(xk+ei,)+args) - apply(f,(xk-ei,)+args))))/(2.0*epsilon)
        ei[k] = 0.0
    return jac

def check_grad(func, grad, x0, *args):
    approx_grad = approx_fprime(x0,func,_epsilon,*args)
    print "Finite difference gradient ", approx_grad
    analytic_grad = grad(x0,*args)
    print "Analytic gradient ", analytic_grad
    differencenorm = sqrt(sum(approx_grad-analytic_grad)**2)
    print "Norm of difference is ", differencenorm
    return differencenorm 

def approx_fhess_p(x0,p,fprime,epsilon,*args):
    f2 = apply(fprime,(x0+epsilon*p,)+args)
    f1 = apply(fprime,(x0,)+args)
    return (f2 - f1)/epsilon

def safe_res(f,x,args):
    """
    Applies f to x.
    Returns f(x) and cost = sum(f(x)**2).
    
    In the case that cost = NaN, returns cost = inf.
    In the case of an exception, returns res = None, cost = inf.
    """
    try:
        res = asarray(apply(f,(x,)+args))
        cost = sum(res**2)
    except (SloppyCell.Utility.SloppyCellException,OverflowError):
        res = None
        cost = scipy.inf
    if scipy.isnan(cost): cost = scipy.inf
    return res, cost

def safe_fprime(fprime,x,args):
    """
    Applies fprime to x.  
    Returns j and exit code.  For nonzero exit codes, j is returned as None.
    
    Exit code 0: No errors.
    Exit code 3: Jacobian contains NaN or inf.
    Exit code 4: Exception in Jacobian calculation.
    """
    try:
        j = asarray(apply(fprime,(x,)+args))
        err = 0
    except SloppyCell.Utility.SloppyCellException:
        j = None
        err = 4
    if j is not None:
      if ( scipy.isnan(j).any() or scipy.isinf(j).any() ):
        j = None
        err = 3
    return j, err

def fmin_lm(f, x0, fprime=None, args=(), avegtol=1e-5, epsilon=_epsilon,
              maxiter=None, full_output=0, disp=1, retall=0, lambdainit = None, 
              jinit = None, trustradius = 1.0):
    """Minimizer for a nonlinear least squares problem. Allowed to
    have more residuals than parameters or vice versa.
    f : residual function (function of parameters)
    fprime : derivative of residual function with respect to parameters.
             Should return a matrix (J) with dimensions  number of residuals
             by number of parameters.
    x0 : initial parameter set
    avegtol : convergence tolerance on the gradient vector
    epsilon : size of steps to use for finite differencing of f (if fprime
              not passed in)
    maxiter : maximum number of iterations 
    full_output : 0 to get only the minimum set of parameters back
                  1 if you also want the best parameter set, the 
                  lowest value of f, the number of function calls, 
                  the number of gradient calls, the convergence flag,
                  the last Marquardt parameter used (lambda), and the 
                  last evaluation of fprime (J matrix)
    disp : 0 for no display, 1 to give cost at each iteration and convergence
           conditions at the end
    retall : 0 for nothing extra to be returned, 1 for all the parameter 
             sets during the optimization to be returned 
    lambdainit : initial value of the Marquardt parameter to use (useful if
                 continuing from an old optimization run
    jinit : initial evaluation of the residual sensitivity matrix (J).
    trustradius : set this to the maximum move you want to allow in a single
                  parameter direction. 
                  If you are using log parameters, then setting this
                  to 1.0, for example, corresponds to a multiplicative 
                  change of exp(1) = 2.718
    """

    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    xcopy = copy.copy(x0)
    if isinstance(x0,KeyedList) :
        x0 = asarray(x0.values())
    else :
        x0 = asarray(x0)

    if lambdainit != None :
        Lambda = lambdainit
    else :
        Lambda = 1.0e-2
    Mult = 10.0
    n = len(x0)
    func_calls = 0
    grad_calls = 0
    res,currentcost = safe_res(f,x0,args)
    func_calls+=1
    m = res.shape[0]
    if maxiter is None :
        maxiter = 200*n
    niters = 0
    x = x0

    gtol = n*avegtol

    if retall:
            allvecs = [x]
    x1 = x0
    x2 = x0
    d = zeros(n,scipy.float_)
    move = zeros(n,scipy.float_)
    finish = 0
    if jinit!=None :
        j = jinit
    else :
        if app_fprime :
            j = asarray(apply(approx_fprime2,(x,f,epsilon)+args))
            func_calls = func_calls + 2*len(x)
        else :
            j,err = safe_fprime(fprime,x,args)
            if err:
                finish = err
            grad_calls+=1

    # NOTE: Below is actually *half* the gradient (because
    # we define the cost as the sum of squares of residuals)
    # However the equations defining the optimization move, dp, 
    # are  2.0*J^tJ dp = -2.0*J^t r, where r is the residual
    # vector; therefore, the twos cancel both sides 
    if j is not None: grad = mat(res)*mat(j)

    while (niters<maxiter) and (finish == 0):
    # note: grad, res and j will be available from the end of the
    # last iteration. They just need to be computed the zeroth
    # time aswell (above)

        lmh = mat(transpose(j))*mat(j)
        # use more accurate way to get e-vals/dirns
        #[u,s,v] = scipy.linalg.svd(lmh)
        [u,ssqrt,vt] = scipy.linalg.svd(j)
        # want n singular values even if m<n and we have
        # more parameters than data points.
        if (len(ssqrt) == n) :
            s = ssqrt**2
        elif (len(ssqrt)<n) :
            s = zeros((n,),scipy.float_)
            s[0:len(ssqrt)] = ssqrt**2
        #print "s is (in original) ", s
        #rhsvect = -mat(transpose(u))*mat(transpose(grad))

        rhsvect = -mat(vt)*mat(transpose(grad))
        rhsvect = asarray(rhsvect)[:,0]
        move = abs(rhsvect)/(s+Lambda*scipy.ones(n)+1.0e-30*scipy.ones(n))
        move = list(move)
        maxindex = move.index(max(move))
        move = asarray(move)

        if max(move) > trustradius :
            Lambda = Mult*(1.0/trustradius*abs(rhsvect[maxindex])-s[maxindex])
            #print " Increasing lambda to ", Lambda
        # now do the matrix inversion

        for i in range(0,n) :
            if (s[i]+Lambda) < 1.0e-30 :
                d[i] = 0.0
            else :
                d[i] = 1.0/(s[i]+Lambda)
            move[i] = d[i]*rhsvect[i]
        move = asarray(move)
        # move = asarray(mat(transpose(v))*mat(transpose(mat(move))))[:,0]
        move = asarray(mat(transpose(vt))*mat(transpose(mat(move))))[:,0]
        # print move
        x1 = x + move
        moveold = move[:]
        
        for i in range(0,n) :
            if (s[i]+Lambda/Mult) < 1.0e-30 :
                d[i] = 0.0
            else :
                d[i] = 1.0/(s[i]+Lambda/Mult)
            move[i] = d[i]*rhsvect[i]
        move = asarray(mat(transpose(vt))*mat(transpose(mat(move))))[:,0]

        x2 = x + asarray(move)
        _,currentcost = safe_res(f,x,args)
        func_calls+=1
        res2,costlambdasmaller = safe_res(f,x2,args)
        func_calls+=1
        res1,costlambda = safe_res(f,x1,args)
        func_calls+=1
        if disp :
            print 'Iteration number', niters
            print 'Current cost', currentcost
            print "Move 1 gives cost of" , costlambda
            print "Move 2 gives cost of ", costlambdasmaller
            #fp = open('LMoutfile','a')
            #fp.write('Iteration number ' + niters.__str__() + '\n')
            #fp.write('Current cost ' + currentcost.__str__() + '\n')
            #fp.write('Move 1 gives cost of ' + costlambda.__str__() + '\n')
            #fp.write('Move 2 gives cost of ' + costlambdasmaller.__str__() + '\n')
            #fp.close()

        oldcost = currentcost
        oldres = res
        oldjac = j

        if costlambdasmaller <= currentcost :
            xprev = x[:]
            Lambda = Lambda/Mult
            x = x2[:]
            if retall:
                allvecs.append(x)
            currentcost = costlambdasmaller
            if app_fprime :
                j = asarray(apply(approx_fprime2,(x2,f,epsilon)+args))
                func_calls = func_calls + 2*len(x2)
            else :
                j,err = safe_fprime(fprime,x2,args)
                if err:
                    x = xprev[:]
                    finish = err
                grad_calls+=1
            if j is not None: grad = mat(res2)*mat(j)
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
        elif costlambda <= currentcost :
            xprev = x[:]
            currentcost = costlambda
            x = x1[:]
            move = moveold[:]
            if retall:
                allvecs.append(x)
            if app_fprime :
                j = asarray(apply(approx_fprime2,(x1,f,epsilon)+args))
                func_calls = func_calls + 2*len(x1)
            else :
                j,err = safe_fprime(fprime,x1,args)
                if err:
                    x = xprev[:]
                    finish = err
                grad_calls+=1
            if j is not None: grad = mat(res1)*mat(j)
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
        else :
            Lambdamult = Lambda
            costmult = costlambda
            piOverFour = .78539816339744825
            NTrials = 0
            NTrials2 = 0
            move = moveold[:]
            while (costmult > currentcost) and (NTrials < 10) :
                num = -scipy.dot(grad,move)[0]
                den = scipy.linalg.norm(grad)*scipy.linalg.norm(move)
                gamma = scipy.arccos(num/den)
                NTrials = NTrials+1
                # was (gamma>piOverFour) below but that doens't
                # make much sense to me. I don't think you should
                # cut back on a given step, I think the trust
                # region strategy is more successful
                if (gamma > 0) :
                    Lambdamult = Lambdamult*Mult
                    for i in range(0,n) :
                        if s[i]+Lambdamult < 1.0e-30 :
                            d[i] = 0.0
                        else :
                            d[i] = 1.0/(s[i]+Lambdamult)
                        move[i] = d[i]*rhsvect[i]
                    move = asarray(mat(transpose(vt))*mat(transpose(mat(move))))[:,0]
                    x1 = x + move
                    res1,costmult = safe_res(f,x1,args)
                    func_calls+=1
                    
                else :
                    NTrials2 = 0
                    while (costmult > currentcost) and (NTrials2 < 10) :
                        NTrials2 = NTrials2 + 1
                        if disp == 1:
                            print " Decreasing stepsize "
                        move = (.5)**NTrials2*moveold
                        x1 = x + asarray(move)
                        res1,costmult = safe_res(f,x1,args)
                        func_calls+=1

            if (NTrials==10) or (NTrials2==10) :
                if disp == 1:
                    print " Failed to converge"
                finish = 1
            else :
                xprev = x[:]
                x = x1[:]
                if retall:
                    allvecs.append(x)
                Lambda = Lambdamult
                if app_fprime :
                    j = asarray(apply(approx_fprime2,(x,f,epsilon)+args))
                    func_calls = func_calls + 2*len(x)
                else :
                    j,err = safe_fprime(fprime,x,args)
                    if err:
                        x = xprev[:]
                        finish = err
                    grad_calls+=1

                if j is not None: grad = mat(res1)*mat(j)
                currentcost = costmult
                if sum(abs(2.0*grad), axis=None) < gtol :
                    finish = 2
        niters = niters + 1

        # see if we need to reduce the trust region
        newmodelval = oldres+asarray(mat(oldjac)*mat(transpose(mat(move))))[:,0]
        oldmodelval = oldres
        #print oldcost-sum(newmodelval**2)
        #print trustradius
        if ((oldcost-sum(newmodelval**2))>1.0e-16) :
            ratio = (oldcost-currentcost)/(oldcost-sum(newmodelval**2))
            if ratio < .25 :
                trustradius = trustradius/2.0
            if ratio >.25 and ratio<=.75 :
                trustradius = trustradius
            if ratio > .75 and trustradius<10.0 :
                trustradius = 2.0*trustradius

        #save(x,'currentParamsLM')

    if disp :
        if (niters>=maxiter) and (finish != 2) :
            print " Current function value: %f" % currentcost
            print " Iterations: %d" % niters
            print " Function evaluations: %d" % func_calls
            print " Gradient evaluations: %d" % grad_calls
            print " Maximum number of iterations exceeded with no convergence "
        if (finish == 2) :
            print " Optimization terminated successfully."
            print " Current function value: %f" % currentcost
            print " Iterations: %d" % niters
            print " Function evaluations: %d" % func_calls
            print " Gradient evaluations: %d" % grad_calls
        if (finish == 3) :
            print " Optimization aborted: Jacobian contains nan or inf."
            print " Current function value: %f" % currentcost
            print " Iterations: %d" % niters
            print " Function evaluations: %d" % func_calls
            print " Gradient evaluations: %d" % grad_calls
        if (finish == 4) :
            print " Optimization aborted: Exception in Jacobian calculation."
            print " Current function value: %f" % currentcost
            print " Iterations: %d" % niters
            print " Function evaluations: %d" % func_calls
            print " Gradient evaluations: %d" % grad_calls

    if isinstance(xcopy,KeyedList) :
        xcopy.update(x)
    else :
        xcopy = x

    if full_output:
        retlist = xcopy, currentcost, func_calls, grad_calls, finish, Lambda, j
        if retall:
            retlist += (allvecs,)
    else :
        retlist = xcopy
        if retall :
            retlist = (xcopy,allvecs)

    return retlist

def fmin_lmNoJ(fcost, x0, fjtj, args=(), avegtol=1e-5, epsilon=_epsilon,
              maxiter=None, full_output=0, disp=1, retall=0, trustradius=1.0):
    """Minimizer for a nonlinear least squares problem. Allowed to
    have more residuals than parameters or vice versa
    fcost : the cost function (*not* the residual function)
    fjtj : this function must return back an ordered pair, the first entry 
           is the gradient of the cost and the second entry is the Levenberg
           Marquardt (LM) approximation to the cost function. 
           NOTE: If the cost function = 1/2 * sum(residuals**2) then 
           the LM approximation is the matrix matrix product J^t J 
           where J = derivative of residual function with respect to parameters. 
           However if cost = k*sum(residuals**2) for some constant k, then 
           the LM approximation is 2*k*J^t J, so beware of this factor!!!
    x0 : initial parameter set
    avegtol : convergence tolerance on the gradient vector
    epsilon : size of steps to use for finite differencing of f (if fprime
              not passed in)
    maxiter : maximum number of iterations 
    full_output : 0 to get only the minimum set of parameters back
                  1 if you also want the best parameter set, the 
                  lowest value of f, the number of function calls, 
                  the number of gradient calls, the convergence flag,
                  the last Marquardt parameter used (lambda), and the 
                  last evaluation of fprime (J matrix)
    disp : 0 for no display, 1 to give cost at each iteration and convergence
           conditions at the end
    retall : 0 for nothing extra to be returned, 1 for all the parameter 
             sets during the optimization to be returned 
    trustradius : set this to the maximum move you want to allow in a single
                  parameter direction. 
                  If you are using log parameters, then setting this
                  to 1.0, for example, corresponds to a multiplicative 
                  change of exp(1) = 2.718
    
    
    This version requires fjtj to pass back an ordered pair with 
    a gradient evaluation of the cost and JtJ,  but not a function for J. 
    This is important in problems when there is many residuals and J is too 
    cumbersome to compute and pass around, but JtJ is a lot "slimmer".  """
    
    xcopy = copy.copy(x0)
    if isinstance(x0,KeyedList) :
        x0 = asarray(x0.values())
    else :
        x0 = asarray(x0)
    Lambda = 1.0e-02
    Mult = 10.0
    n = len(x0)
    func_calls = 0
    grad_calls = 0
    if maxiter==None :
        maxiter = 200*n
    niters = 0
    x = x0

    gtol = n*avegtol

    if retall:
        allvecs = [x]
    x1 = x0
    x2 = x0
    d = zeros(n,scipy.float_)
    move = zeros(n,scipy.float_)
    finish = 0

    grad, lmh = apply(fjtj,(x,))
    grad_calls+=1


    while (niters<maxiter) and (finish == 0):
        # estimate what Lambda should be

        [u,s,v] = scipy.linalg.svd(lmh)
        #print "s is (in NoJ) ", s
        #s,u = scipy.linalg.eig(lmh)
        #s = real(s)
        #u = real(u)
        oldlmh = lmh[:,:]
        oldgrad = grad[:]

        rhsvect = -scipy.dot(transpose(u),grad)
        # rhsvect = asarray(rhsvect)[:,0]
        move = abs(rhsvect)/(s+Lambda*ones(n)+1.0e-30*ones(n))
        move = list(move)
        maxindex = move.index(max(move))
        move = asarray(move)
        if max(move) > trustradius :
            Lambda = Mult*(1.0/trustradius*abs(rhsvect[maxindex])-s[maxindex])
            #print " Increasing lambda to ", Lambda

        ## lmhreg = lmh + Lambda*eye(n,n,typecode=scipy.float_)
        ## [u,s,v] = scipy.linalg.svd(lmhreg)
        rhsvect = -scipy.dot(transpose(u),grad)
        # rhsvect = asarray(rhsvect)[:,0]

        for i in range(0,len(s)) :
            if (s[i]+Lambda) < 1.0e-30 :
                d[i] = 0.0
            else :
                d[i] = 1.0/(s[i]+Lambda)
            move[i] = d[i]*rhsvect[i]
        move = asarray(move)
        move = dot(asarray(u),move)
        x1 = x + move
        moveold = move[:]


        for i in range(0,len(s)) :
            if (s[i]+Lambda/Mult) < 1.0e-30 :
                d[i] = 0.0
            else :
                d[i] = 1.0/(s[i]+Lambda/Mult)
            move[i] = d[i]*rhsvect[i]
        move = asarray(move)
        move = dot(asarray(u),move)
        x2 = x + asarray(move)

        currentcost = apply(fcost,(x,))
        oldcost = currentcost
        func_calls+=1

        try:
            costlambdasmaller = apply(fcost,(x2,))
        except SloppyCell.Utility.SloppyCellException:
            costlambdasmaller = scipy.inf
        func_calls+=1

        try:
            costlambda = apply(fcost,(x1,))
        except SloppyCell.Utility.SloppyCellException:
            costlambda = scipy.inf
        func_calls+=1
        if disp :
            print 'Iteration number', niters
            print 'Current cost', currentcost
            print "Move 1 gives cost of" , costlambda
            print "Move 2 gives cost of ", costlambdasmaller

            #fp = open('LMoutfile','a')
            #fp.write('Iteration number ' + niters.__str__() + '\n')
            #fp.write('Current cost ' + currentcost.__str__() + '\n')
            #fp.write('Move 1 gives cost of ' + costlambda.__str__() + '\n')
            #fp.write('Move 2 gives cost of ' + costlambdasmaller.__str__() + '\n')
            #fp.close()

        if costlambdasmaller <= currentcost :
            Lambda = Lambda/Mult
            x = x2[:]
            if retall:
                allvecs.append(x)
            currentcost = costlambdasmaller
            grad, lmh = apply(fjtj,(x2,))
            grad_calls+=1

            #if scipy.linalg.norm(asarray(grad)) < avegtol :
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
        elif costlambda <= currentcost :
            currentcost = costlambda
            move = moveold[:]
            x = x1[:]
            if retall:
                allvecs.append(x)
            grad, lmh = apply(fjtj,(x1,))
            grad_calls+=1

            # if scipy.linalg.norm(asarray(grad)) < avegtol :
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
        else :
            Lambdamult = Lambda
            costmult = costlambda
            piOverFour = .78539816339744825
            NTrials2 = 0
            NTrials = 0

            while (costmult > currentcost) and (NTrials < 10) :
                # num = -dot(transpose(asarray(grad)),asarray(moveold) )
                # den = scipy.linalg.norm(grad)*scipy.linalg.norm(moveold)
                gamma = .1 # scipy.arccos(num/den)
                NTrials = NTrials+1
                if (gamma > 0) :
                    Lambdamult = Lambdamult*Mult

                    for i in range(0,len(s)) :
                        if s[i] + Lambdamult < 1.0e-30 :
                            d[i] = 0.0
                        else :
                            d[i] = 1.0/(s[i] + Lambdamult)
                        move[i] = d[i]*rhsvect[i]
                    move = asarray(move)
                    move = dot(asarray(u),move)
                    x1 = x + asarray(move)

                    func_calls+=1
                    costmult = apply(fcost,(x1,))
                else :
                    NTrials2 = 0
                    while (costmult > currentcost) and (NTrials2 < 10) :
                        NTrials2 = NTrials2 + 1
                        if disp :
                            print " Decreasing stepsize "
                        move = (.5)**NTrials2*moveold
                        x1 = x + asarray(moveold)
                        func_calls+=1
                        costmult = apply(fcost,(x1,))

            if (NTrials==10) or (NTrials2==10) :
                if disp :
                    print " Failed to converge"
                finish = 1
            else :
                x = x1[:]
                if retall:
                    allvecs.append(x)
                Lambda = Lambdamult
                grad, lmh = apply(fjtj,(x1,))
                grad_calls+=1
                currentcost = costmult
                # if scipy.linalg.norm(grad) < avegtol :
                if sum(abs(2.0*grad), axis=None) < gtol :
                    finish = 2
        niters = niters + 1

        # see if we need to reduce the trust region, compare the actual change in
        # cost to the linear and quadratic change in cost
        model_change = scipy.dot(scipy.transpose(oldgrad),move) + \
        .5*scipy.dot(scipy.transpose(move),scipy.dot(oldlmh,move) )
        #print oldcost-sum(newmodelval**2)
        #print trustradius
        if model_change>1.0e-16 :
            ratio = (oldcost-currentcost)/(model_change)
            if ratio < .25 :
                trustradius = trustradius/2.0
            if ratio >.25 and ratio<=.75 :
                trustradius = trustradius
            if ratio > .75 and trustradius<10.0 :
                trustradius = 2.0*trustradius

        #save(x,'currentParamsLM')

    if disp :
        if (niters>=maxiter) and (finish != 2) :
            print "         Current function value: %f" % currentcost
            print "         Iterations: %d" % niters
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
            print " Maximum number of iterations exceeded with no convergence "
        if (finish == 2) :
            print "Optimization terminated successfully."
            print "         Current function value: %f" % currentcost
            print "         Iterations: %d" % niters
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls

    if isinstance(xcopy,KeyedList) :
        xcopy.update(x)
    else :
        xcopy = x

    if full_output:
        retlist = xcopy, currentcost, func_calls, grad_calls, finish, Lambda, lmh
        if retall:
            retlist += (allvecs,)
    else:
        retlist = xcopy
        if retall:
            retlist = (xcopy, allvecs)

    return retlist

def solve_lmsys(Lambda,s,g,rhsvect,currentcost,n) :
    d = zeros(n,scipy.float_)
    move = zeros(n,scipy.float_)
    for i in range(0,n) :
        if s[i] < 1.0e-20 :
            d[i] = 0.0
        else :
            d[i] = 1.0/(s[i])
        move[i] = d[i]*rhsvect[i]
    return move

def fmin_lm_scale(f, x0, fprime=None, args=(), avegtol=1e-5, epsilon=_epsilon,
              maxiter=None, full_output=0, disp=1, retall=0,trustradius=1.0):
    """
    Minimizer for a nonlinear least squares problem. Allowed to
    have more residuals than parameters or vice versa. 

    f : residual function (function of parameters)
    fprime : derivative of residual function with respect to parameters.
             Should return a matrix (J) with dimensions  number of residuals
             by number of parameters.
    x0 : initial parameter set
    avegtol : convergence tolerance on the gradient vector
    epsilon : size of steps to use for finite differencing of f (if fprime
              not passed in)
    maxiter : maximum number of iterations 
    full_output : 0 to get only the minimum set of parameters back
                  1 if you also want the best parameter set, the 
                  lowest value of f, the number of function calls, 
                  the number of gradient calls, the convergence flag,
                  the last Marquardt parameter used (lambda), and the 
                  last evaluation of fprime (J matrix)
    disp : 0 for no display, 1 to give cost at each iteration and convergence
           conditions at the end
    retall : 0 for nothing extra to be returned, 1 for all the parameter 
             sets during the optimization to be returned 
    trustradius : set this to the maximum length of move you want. 
                  If you are using log parameters, then setting this
                  to 1.0, for example, corresponds to a multiplicative 
                  change of exp(1) = 2.718 if the move is along a single
                  parameter direction
    
    This version is scale invariant. This means that under a change of 
    scale of the parameters the direction the optimizer chooses to move
    in does not change. To achieve this, we don't use a Marquardt
    parameter to impose a trust region but rather take the infinite trust 
    region step and just cut it back to the length given in the variable
    trustradius. """
    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    xcopy = copy.copy(x0)
    if isinstance(x0,KeyedList) :
        x0 = asarray(x0.values())
    else :
        x0 = asarray(x0)

    Lambda = 1.0e-02
    Mult = 10.0
    n = len(x0)
    func_calls = 0
    grad_calls = 0
    res = asarray(apply(f,(x0,)))
    m = res.shape[0]
    if maxiter is None :
        maxiter = 200*n
    niters = 0
    x = x0

    gtol = n*avegtol

    if retall:
            allvecs = [x]
    x1 = x0
    x2 = x0
    d = zeros(n,scipy.float_)
    move = zeros(n,scipy.float_)
    finish = 0

    if app_fprime :
        j = asarray(apply(approx_fprime2,(x,f,epsilon)+args))
        func_calls = func_calls + 2*len(x)
    else :
        j = asarray(apply(fprime,(x,)))
        grad_calls+=1
    res = asarray(apply(f,(x,)))
    func_calls+=1
    grad = mat(res)*mat(j)

    while (niters<maxiter) and (finish == 0):
    # note: grad, res and j will be available from the end of the
    # last iteration. They just need to be computed the zeroth
    # time aswell (above)

        lmh = mat(transpose(j))*mat(j)
        # use more accurate way to get e-vals/dirns
        #[u,s,v] = scipy.linalg.svd(lmh)
        [u,ssqrt,vt] = scipy.linalg.svd(j)
        # want n singular values even if m<n and we have
        # more parameters than data points.
        if (len(ssqrt) == n) :
            s = ssqrt**2
        elif (len(ssqrt)<n) :
            s = zeros((n,),scipy.float_)
            s[0:len(ssqrt)] = ssqrt**2
        #rhsvect = -mat(transpose(u))*mat(transpose(grad))
        rhsvect = -mat(vt)*mat(transpose(grad))
        rhsvect = asarray(rhsvect)[:,0]


        currentcost = sum(asarray(apply(f,(x,)))**2)
        g = asarray(grad)[0,:]
        Lambda = 0
        move = solve_lmsys(Lambda,s,g,rhsvect,currentcost,n)

        move = asarray(move)
        move = asarray(mat(transpose(vt))*mat(transpose(mat(move))))[:,0]
        unitmove = move/(scipy.linalg.norm(move))
        move1 = unitmove*trustradius

        # print move
        x1 = x + move1
        
        move2 = unitmove*trustradius*Mult
        x2 = x + asarray(move2)

        func_calls+=1
        try:
            res2 = asarray(apply(f,(x2,)))
            costlambdasmaller = sum(res2**2)
        except SloppyCell.Utility.SloppyCellException:
            costlambdasmaller = scipy.inf
        func_calls+=1
        try:
            res1 = asarray(apply(f,(x1,)))
            costlambda = sum(res1**2)
        except SloppyCell.Utility.SloppyCellException:
            costlambda = scipy.inf
        func_calls+=1
        if disp :    
            print "Cost is ", currentcost
            print "Iteration is", niters

        oldcost = currentcost
        oldres = res
        oldjac = j

        if costlambdasmaller <= currentcost :
            trustradius = trustradius*Mult
            x = x2
            if retall:
                allvecs.append(x)
            currentcost = costlambdasmaller
            if app_fprime :
                j = asarray(apply(approx_fprime2,(x2,f,epsilon)+args))
                func_calls = func_calls + 2*len(x2)
            else :
                j = asarray(apply(fprime,(x2,)))
                grad_calls+=1
            grad = mat(res2)*mat(j)
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
            move = move2
        elif costlambda <= currentcost :
            currentcost = costlambda
            x = x1
            if retall:
                allvecs.append(x)
            if app_fprime :
                j = asarray(apply(approx_fprime2,(x1,f,epsilon)+args))
                func_calls = func_calls + 2*len(x1)
            else :
                j = asarray(apply(fprime,(x1,)))
                grad_calls+=1

            grad = mat(res1)*mat(j)
            if sum(abs(2.0*grad), axis=None) < gtol :
                finish = 2
            move = move1
        else :
            trustradmult = trustradius
            costmult = costlambda
            NTrials = 0
            move = unitmove
            while (costmult > currentcost) and (NTrials < 100) :
                while (costmult > currentcost) and (NTrials < 100) :
                    NTrials = NTrials + 1
                    #print " Decreasing stepsize "
                    trustradmult = trustradmult/2.0
                    move = move*trustradmult
                    x1 = x + asarray(move)
                    res1 = asarray(apply(f,(x1,)))
                    func_calls+=1
                    costmult = sum(res1**2)

            if (NTrials==100) :
                if disp :
                    print " Failed to converge"
                finish = 1
            else :
                x = x1
                if retall:
                    allvecs.append(x)
                trustradius = trustradmult
                if app_fprime :
                    j = asarray(apply(approx_fprime2,(x,f,epsilon)+args))
                    func_calls = func_calls + 2*len(x)
                else :
                    j = asarray(apply(fprime,(x,)))
                    grad_calls+=1

                grad = mat(res1)*mat(j)
                currentcost = costmult
                if sum(abs(2.0*grad), axis=None) < gtol :
                    finish = 2
        niters = niters + 1

        # see if we need to reduce the trust region
        newmodelval = oldres+asarray(mat(oldjac)*mat(transpose(mat(move))))[:,0]
        oldmodelval = oldres
        #print oldcost-sum(newmodelval**2)
        #print trustradius
        if ((oldcost-sum(newmodelval**2))>1.0e-16) :
            ratio = (oldcost-currentcost)/(oldcost-sum(newmodelval**2))
            if ratio < .25 :
                trustradius = trustradius/2.0
            if ratio >.25 and ratio<=.75 :
                trustradius = trustradius
            if ratio > .75 and trustradius<10.0 :
                trustradius = 2.0*trustradius

    if disp :
        if (niters>=maxiter) and (finish != 2) :
            print "         Current function value: %f" % currentcost
            print "         Iterations: %d" % niters
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
            print " Maximum number of iterations exceeded with no convergence "
        if (finish == 2) :
            print "Optimization terminated successfully."
            print "         Current function value: %f" % currentcost
            print "         Iterations: %d" % niters
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
            
    if isinstance(xcopy,KeyedList) :
        xcopy.update(x)
    else :
        xcopy = x

    if full_output:
        retlist = xcopy, currentcost, func_calls, grad_calls, finish, Lambda, j 
        if retall:
            retlist += (allvecs,)
    else:
        retlist = xcopy
        if retall:
            retlist = (xcopy, allvecs)

    return retlist
