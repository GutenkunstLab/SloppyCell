import scipy
import time
try:
    import SloppyCell.Plotting as Plotting
except ImportError:
    pass

def C(p, origMatrix, weightsRDP=0.,weights2D=0.,weightsLS=0.,weightPR=0.,weightPriors=0.,*args,**kwargs):
    """
    p is list {a_{ij}, i: 1..n, j: i+1..n}
    n is size of total matrix
    """
    n = origMatrix.shape[0]
    OrthMatrix = ProcessHalfMatrix(p)
    newMat = transformMatrix(origMatrix,OrthMatrix)
    cost=0.
    if weightsRDP > 0.:
        cost += weightsRDP*sumRowDotProdsOLD(newMat)
    if weights2D > 0.:
        cost += weights2D*sum2Determinants(newMat)
    if weightsLS > 0.:
        cost += weightsLS*sumLogSpacings(newMat)
    if weightPR>0.:
        cost += weightPR/(ParticipationRatio(OrthMatrix)**2.)
##         cost += weightPR*(1.-ParticipationRatio(OrthMatrix)/n)
    if weightPriors>0.:
        cost +=weightPR*calcPriors(p)
    return cost

def calcPriors(paramsList, pOpt=0.,pSigma=10.):
    priorCost=0.
    for param in paramsList:
        priorCost += ((param-pOpt)/pSigma)**2.
    return priorCost

def sumRowDotProdsOLD(origMatrix):
    """
    Makes more sense to use on Jacobian than on Hessian.
    """
    rowNormalized = normRows(origMatrix)
    n = origMatrix.shape[0]
##    sumDotProds = sum([abs((scipy.dot(rowNormalized[i],rowNormalized[i+1])))**(1./2.) for i in range(n-1)])
    sumDotProds = sum([1.-(scipy.dot(rowNormalized[i],rowNormalized[i+1]))**2. for i in range(n-1)])
    return sumDotProds

def sumRowDotProdsNEW(origMatrix):
    """
    Makes more sense to use on Jacobian than on Hessian.
    """
    rowNormalized = normRows(origMatrix)
    n = rowNormalized.shape[0]
    sumDotProds = sumAllDotProds(rowNormalized[0:n/2]) + sumAllDotProds(rowNormalized[n/2:n])
    return sumDotProds

def sumAllDotProds(origMatrix):
    n = origMatrix.shape[0]
    sumDotProds=0.
    for ii in range(n):
        for jj in range(ii+1,n):
            sumDotProds += 1.-scipy.dot(origMatrix[ii],origMatrix[jj])**2.
    return sumDotProds

def sum2Determinants(matrixToSum):
    n = matrixToSum.shape[0]/2
    det1 = scipy.linalg.det(matrixToSum[0:n,0:n])
    det2 = scipy.linalg.det(matrixToSum[n:n*2,n:n*2])
    return scipy.absolute(det1)+scipy.absolute(det2)

def sumLogSpacings(matrixToSum):
    n = matrixToSum.shape[0]/2
    sv1 = scipy.linalg.svdvals(matrixToSum[0:n,0:n])
    sv2 = scipy.linalg.svdvals(matrixToSum[n:n*2,n:n*2])
    return sum(sv1[1:n]/sv1[0:n-1])+sum(sv2[1:n]/sv2[0:n-1])

def ParticipationRatio(origMatrix):
    return scipy.sum(scipy.sum(origMatrix**4))

def ProcessFullMatrix(p):
    """
    this is used if parameters define all elements of a matrix
    """
    n = scipy.sqrt(len(p)).__int__()
    pMat = scipy.reshape(p,(n,n))
    orthMat = scipy.linalg.orth(pMat)
    return orthMat

def ProcessHalfMatrix(p):
    """
    this is used if parameters define upper right triangular portion of
    skew-symmetric matrix.
    Then performs a Cayley transformation of this skew-symmetric matrix.
    """
    l = scipy.size(p)
    n = scipy.round(((1.+scipy.sqrt(1+8.*l))/2.)).__int__()
    Mfull = scipy.zeros((n,n),'d')
    placeCounter=0
    for i in range(n-1):
        Mfull[i,(i+1):n] = p[placeCounter:(placeCounter+n-i-1)]
        placeCounter += n-i-1
    Mfull = Mfull - scipy.transpose(Mfull)
    IMat = scipy.eye(n,'d')
    return scipy.dot((IMat-Mfull),scipy.linalg.inv(IMat+Mfull))
##     return scipy.linalg.expm(Mfull)

def ProcessQuarterMatrix(p):
    """
    this is used if parameters define just upper right and lower left
    quadrants of antisymmetric matrix.
    """
    n= scipy.sqrt(scipy.size(p))
    Mu = scipy.resize(p, (n,n))
    Mfull = scipy.bmat([[scipy.zeros((n,n)),Mu],[-scipy.transpose(Mu),scipy.zeros((n,n))]])
    return scipy.linalg.expm(Mfull)

def BestMatrix(origMatrix, p=None, weightsRDP=1., weights2D=0., weightsLS=0., weightPR=0., weightPriors=0., seed=None, *args, **kwargs):
    if seed is not None:
        scipy.random.seed(seed)
    n = origMatrix.shape[0]
    if p is None:
        p = (scipy.random.random(n*(n-1)/2)-0.5)
    pOpt = scipy.optimize.fmin(C, p, args=(origMatrix,weightsRDP, weights2D, weightsLS, weightPR, weightPriors), *args, **kwargs)
    return C(pOpt, origMatrix, weightsRDP, weights2D, weightsLS, weightPR, weightPriors), pOpt

def ManyBest(origMatrix, numberTries=10, p=None,numOpts=10,weightsRDP=1., weights2D=0., weightsLS=0., weightPR=0., weightPriors=0., seed=None, *args, **kwargs):
    return [OptimizeManyTimes(origMatrix,p,numOpts,weightsRDP, weights2D, weightsLS, weightPR, weightPriors, seed, *args, **kwargs) for i in range(numberTries)]
    
def OptimizeManyTimes(origMatrix,p=None,numOpts=10,weightsRDP=1., weights2D=0., weightsLS=0., weightPR=0., weightPriors=0., seed=None, *args, **kwargs):
    C,pOpt = BestMatrix(origMatrix,p,weightsRDP,weights2D,weightsLS,weightPR,weightPriors,seed,*args,**kwargs)
    for i in range(numOpts-1):
        C,pOpt = BestMatrix(origMatrix,pOpt,weightsRDP,weights2D,weightsLS,weightPR,weightPriors,seed,*args,**kwargs)
    return C, pOpt

def getGammas(numExps,distWidth=1.):
    epsilons = scipy.random.random(numExps)-0.5
    gammas = 1.+epsilons
    return gammas

def getAmounts(numExps,distWidth=1.):
#    amounts=scipy.ones(numExps,'d')
    amounts = scipy.random.random(numExps)-0.5
    amounts = 1.+amounts
    return amounts

def netRadiation(gammas,amounts,time):
    netR = amounts*scipy.exp(-gammas*time)
    return scipy.sum(netR)

def getHessFull(gammas=None,amounts=None,numExps=3):
    if gammas is None:
        gammas = getGammas(numExps)
    if amounts is None:
        amounts = getAmounts(numExps)

    numExps=scipy.size(gammas)

    hessgg = getHessGG(gammas)

    hessAA = getHessAA(amounts,gammas)

    numerAg1 = -scipy.outer(amounts,amounts*gammas)
    denomAg1 = scipy.outer(scipy.ones(numExps),gammas)
    denomAg2 = scipy.transpose(denomAg1)+denomAg1
    denomAg3 = denomAg2*denomAg2
    hessAg = numerAg1/denomAg3

    numergA1 = -scipy.outer(amounts*gammas,amounts)
    denomgA1 = scipy.outer(gammas,scipy.ones(numExps))
    denomgA2 = scipy.transpose(denomgA1)+denomgA1
    denomgA3 = denomgA2*denomgA2
    hessgA = numergA1/denomgA3

    hessFull = scipy.bmat([[hessAA,hessAg],[hessgA,hessgg]])

    return hessFull

def getHessFullLogTime(gammas=None,amounts=None,numExps=3):
    """use this routine for new paper"""
    if gammas is None:
        gammas = getGammas(numExps)
    if amounts is None:
        amounts = getAmounts(numExps)

    numExps=scipy.size(gammas)

    amountsLess1 = amounts[0:-1]
    gammasLess1 = gammas[0:-1]

    hessgg = getHessGGLogTime(gammas=gammas,amounts=amounts)

    hessAA = getHessAALogTime(amounts=amounts,gammas=gammas)

    numerAg1 = scipy.array(-2.*scipy.outer(amountsLess1,amounts*gammas))
    denomAg1 = scipy.array(scipy.outer(gammasLess1,scipy.ones(numExps)))
    denomAg2 = scipy.array(scipy.outer(scipy.ones(numExps-1),gammas))
    denomAg3 = scipy.array(denomAg1+denomAg2)
    denomAg4 = denomAg2+gammas[-1]
    hessAg = scipy.array(numerAg1/denomAg3)-scipy.array(numerAg1/denomAg4)

    hessgA = scipy.transpose(hessAg)

    hessFull = scipy.array(scipy.bmat([[hessAA,hessAg],[hessgA,hessgg]]))

    return hessFull

def getHessGG(gammas=None,numExps=3):
    if gammas is None:
        numExps=6
        gammas = scipy.array([1.,1.2,1.4,1001.,1005.,997.])

    numExps = scipy.size(gammas)
    numergg1 = 2.*scipy.outer(gammas, gammas)
    denomgg1 = scipy.outer(scipy.ones(numExps),gammas)
    denomgg2 = scipy.transpose(denomgg1)+denomgg1
    denomgg3 = denomgg2*denomgg2*denomgg2
    hessgg = numergg1/denomgg3

    return hessgg

def getHessGGLogTime(gammas=None,amounts=None,numExps=3):
    """use this routine for the new paper"""
    if gammas is None:
        numExps=6
        gammas = scipy.array([1.,1.2,1.4,1001.,1005.,997.])
    if amounts is None:
        amounts = scipy.ones(numExps,'d')

    numExps = scipy.size(gammas)
    numergg1 = scipy.array(2.*scipy.outer(gammas*amounts, gammas*amounts))
    denomgg1 = scipy.array(scipy.outer(scipy.ones(numExps),gammas))
    denomgg2 = scipy.array(scipy.transpose(denomgg1)+denomgg1)
    denomgg3 = scipy.array(denomgg2*denomgg2)
    hessgg = scipy.array(numergg1/denomgg3)

    return hessgg

def getHessAA(amounts=None,gammas=None,numExps=3):
    if amounts is None:
        amounts = getAmounts(numExps)
    if gammas is None:
        gammas = getGammas(numExps)

    numExps = scipy.size(amounts)

    numerAA1 = scipy.outer(amounts,amounts)
    denomAA1 = scipy.outer(scipy.ones(numExps),gammas)
    denomAA2 = scipy.transpose(denomAA1)+denomAA1
    hessAA = numerAA1/denomAA2

    return hessAA

def getHessAALogTime(amounts=None,gammas=None,numExps=3):
    """use this routine for the new paper"""
    if amounts is None:
        amounts = getAmounts(numExps)
    if gammas is None:
        gammas = getGammas(numExps)

    numExps = scipy.size(amounts)
    amountsLess1 = amounts[0:-1]
    gammasLess1 = gammas[0:-1]

    numerAA1 = scipy.array(scipy.outer(amountsLess1,amountsLess1))
    numerAA2 = scipy.array(scipy.outer(gammasLess1+gammas[-1],scipy.ones(numExps-1)))
    numerAA3 = scipy.array(scipy.outer(scipy.ones(numExps-1),gammasLess1+gammas[-1]))
    denomAA1 = scipy.array(scipy.outer(scipy.ones(numExps-1),2*gammas[-1]*gammasLess1))
    denomAA2 = scipy.array(scipy.transpose(denomAA1)+denomAA1)
    hessAA = 2.*scipy.array(numerAA1*scipy.log(numerAA2*numerAA3/denomAA2))

    return hessAA

def getJacobian(times, amounts=None, gammas=None):
    numGammas=0
    numAmounts=0
    if gammas is not None:
        numGammas = scipy.size(gammas)
    if amounts is not None:
        numAmounts = scipy.size(amounts)
    Jacobian = scipy.zeros((numGammas+numAmounts,scipy.size(times)),'d')
    if numAmounts > 0:
        Jacobian[0:numAmounts] = getJacobianA(times,amounts,gammas)
    if numGammas > 0:
        Jacobian[numAmounts:numAmounts+numGammas] = getJacobianG(times,gammas,amounts)
    return Jacobian

def getJacobianLogTime(times, amounts=None, gammas=None):
    numGammas=0
    numAmounts=0
    if gammas is not None:
        numGammas = scipy.size(gammas)
    if amounts is not None:
        numAmounts = scipy.size(amounts)
    JacobianLogTime = scipy.zeros((numGammas+numAmounts,scipy.size(times)),'d')
    if numAmounts > 0:
        JacobianLogTime[0:numAmounts] = getJacobianALogTime(times,amounts,gammas)
    if numGammas > 0:
        JacobianLogTime[numAmounts:numAmounts+numGammas] = getJacobianGLogTime(times,gammas,amounts)
    return JacobianLogTime

def getJacobianLog(times,gammas,amounts=None):
    """This is the routine to use for new paper (with times exponentially
    distributed."""
    numExps = scipy.size(gammas)
    if amounts is None:
        amounts = scipy.ones(numExps,'d')
    jacG = getJacobianLogG(times,gammas,amounts)
    jacA = getJacobianLogA(times,gammas,amounts)
    return scipy.transpose(scipy.array(scipy.bmat([[jacA],[jacG]])))
    

def getJacobianLogG(times,gammas,amounts=None):
    """This is the routine to use for new paper (with times exponentially
    distributed."""
    numExps = scipy.size(gammas)
    if amounts is None:
        amounts = scipy.ones(numExps,'d')
    JacobianLogG = scipy.array([gammas[j]*amounts[j]*(-times)*scipy.exp(-gammas[j]*times) for j in range(numExps)])
    return JacobianLogG

def getJacobianLogA(times,gammas,amounts):
    """This is the routine to use for new paper (with times exponentially
    distributed."""
    numExps = scipy.size(gammas)
    JacobianLogA = scipy.array([amounts[j]*(scipy.exp(-gammas[j]*times)-scipy.exp(-gammas[-1]*times)) for j in range(numExps-1)])
    return JacobianLogA

def getJacobianG(times,gammas,amounts=None):
    numGammas = scipy.size(gammas)
    if amounts is None:
        amounts = scipy.ones(numGammas,'d')
    JacobianG = scipy.array([-amounts[i]*gammas[i]*times*scipy.exp(-gammas[i]*times) for i in range(numGammas)])
    return JacobianG

def getJacobianGLogTime(times,gammas,amounts=None):
    numGammas = scipy.size(gammas)
    if amounts is None:
        amounts = scipy.ones(numGammas,'d')
    JacobianGLogTime = scipy.array([-amounts[i]*gammas[i]*scipy.sqrt(times)*scipy.exp(-gammas[i]*times) for i in range(numGammas)])
    return JacobianGLogTime

def getJacobianA(times,amounts,gammas=None):
    numAmounts = scipy.size(amounts)
    if gammas is None:
        gammas = scipy.zeros(numAmounts,'d')
    JacobianA = scipy.array([amounts[i]*scipy.exp(-gammas[i]*times) for i in range(numAmounts)])
    return JacobianA

def getJacobianALogTime(times,amounts,gammas=None):
    numAmounts = scipy.size(amounts)
    if gammas is None:
        gammas = scipy.zeros(numAmounts,'d')
    JacobianALogTime = scipy.array([(amounts[i]/scipy.sqrt(times))*scipy.exp(-gammas[i]*times) for i in range(numAmounts)])
    return JacobianALogTime

def normRows(matrixToPlot):
    """
    Useful for plotting and otherwise comparing alignment of rows of matrices.
    Be careful that if vec[1] is zero, entire row gets zerod.
    """
    return scipy.array([scipy.sign(vec[1])*vec/scipy.sqrt(scipy.dot(vec,vec)) for vec in matrixToPlot])

def transformMatrix(origMatrix, transformation):
    n, m =origMatrix.shape
    if n == m:
        newMat = scipy.dot(transformation,scipy.dot(origMatrix,scipy.transpose(transformation)))
    else:
        newMat = scipy.dot(transformation,origMatrix)
    return newMat

def findPermutation(permMat):
    """
    In so far as permMat is a permutation matrix, returns the permutation.
    """
    maxs = [(vec.tolist().index(max(vec)), max(vec)) for vec in permMat]
    mins = [(vec.tolist().index(min(vec)), min(vec)) for vec in permMat]
    for ii in range(len(maxs)):
        if maxs[ii][1] < -mins[ii][1]:
            maxs[ii] = mins[ii]
    return maxs

def makePermutationMatrix(permList):
    """
    Takes a list defining the permutation and makes the appropriate matrix.
    """
    permList = scipy.array(permList)
    n = len(permList)
    if 0 not in permList:
        permList = permList - 1
    permMat = scipy.zeros((n,n),'d')
    for ii, jj in enumerate(permList):
        permMat[ii,jj] = 1.
    return permMat

def timeCostEval(mat, p, numEvals=100):
    start = time.time()
    for ii in range(numEvals):
        Cost = C(p,mat)
    stop = time.time()
    print "time per eval:", (stop-start)/numEvals
    return (stop-start)/numEvals

def timeGammas(listNumGammas, numEvals=1000):
    timesToCalc = []
    for numGammas in listNumGammas:
        gammas = getGammas(numGammas)
        jac = getJacobianG(gammas=gammas,times=scipy.arange(500.))
        p = scipy.random.random(numGammas*(numGammas-1)/2)-0.5
        timePerEval = timeCostEval(jac,p,numEvals)
        timesToCalc.append(timePerEval)
    return timesToCalc

def plotLevels(levels,offset=0):
    xs = [offset+0.6,offset+1.4]
    for lvl in levels:
        Plotting.semilogy(xs,[lvl,lvl],'k')
    return

def eVec(k,n):
    eVec = scipy.zeros(n,'d')
    eVec[k] = 1.
    return eVec

def sLSAlongP(origMatrix,params,k,deltaP):
    sLSList = [sumLogSpacings(transformMatrix(origMatrix,ProcessHalfMatrix(params+x*params[k]*eVec(k,len(params))))) for x in scipy.arange(-deltaP,deltaP,deltaP/1000)]
    return sLSList

def PRAlongP(params,k,deltaP):
    n = len(params)
    PRList = [1./ParticipationRatio(ProcessHalfMatrix(params+x*params[k]*eVec(k,n))) for x in scipy.arange(-deltaP,deltaP,deltaP/1000)]
##    PRList = [1./ParticipationRatio(transformMatrix(origMatrix,ProcessHalfMatrix(params+x*params[k]*eVec(k,n)))) for x in scipy.arange(-deltaP,deltaP,deltaP/100)]
##     PRList = [1.-ParticipationRatio(transformMatrix(origMatrix,ProcessHalfMatrix(params+x*params[k]*eVec(k,n))))/n for x in scipy.arange(-deltaP,deltaP,deltaP/100)]
    return PRList

def CostAlongP(origMatrix,params,k,deltaP,weightsRDP=0.,weights2D=0.,weightsLS=0.,weightPR=0.,weightPriors=0.):
    n = len(params)
    CostList = [C(params+x*params[k]*eVec(k,n),origMatrix, weightsRDP,weights2D,weightsLS,weightPR,weightPriors) for x in scipy.arange(-deltaP,deltaP,deltaP/100)]
    return CostList
