import scipy

try:
    import SloppyCell.Plotting as Plotting
except:
    pass

# 1. calculate your jtj matrix
# 2. use calcDmat to get Dmat
# 3. use getLambdaR2 to get lambdaMat, r2Mat
# 4. use calcR2mat_Terms to get ctMat==cross term from r2Mat and
#    dtMat==diagonal term from r2Mat
# 5. use cluster to cluster based on either r2Mat (for unscaled distances) or
#    r2Mat/dtMat (to scale distances by size of individual parameter effects)

def normColumns(mat1):
    n = len(mat1[0])
    mat2 = mat1.copy()
    for ii in range(n):
        mat2[:,ii] = mat2[:,ii]/scipy.sqrt(scipy.dot(mat2[:,ii],mat2[:,ii]))
    return mat2

def colorAlignColumns(mat1):
    n = len(mat1[0])
    mat2 = mat1.copy()
    for ii in range(n-1):
        if scipy.dot(mat2[:,ii],mat2[:,ii+1]) < 0.:
            mat2[:,ii+1] = -1.*mat2[:,ii+1]
    return mat2

def calcCovarianceMat(origMat):
    covarMat = scipy.dot(scipy.transpose(origMat),origMat)
    n = len(covarMat)
    jaa = scipy.outer(scipy.diagonal(covarMat),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(covarMat))
    covarMat = covarMat/scipy.sqrt(jaa)
    covarMat = covarMat/scipy.sqrt(jbb)
    return covarMat

def calcDmat(jtj):
    n = len(jtj)
    Dmat = scipy.zeros((n,n),'d')
    for alpha in range(n):
        for beta in range(alpha,n):
            Dmat[alpha,beta] = (jtj[alpha,alpha]-jtj[beta,beta])**2 + 4.*(jtj[alpha,beta]**2)
    Dmat = Dmat + scipy.transpose(Dmat)
    for alpha in range(n):
        Dmat[alpha,alpha] = 4.*(jtj[alpha,alpha]**2)
    return Dmat

def calcR2mat_Terms(jtj,lambdaMat):
    n = len(jtj)
    jaa = scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj))
    crosstermMat = 2.*lambdaMat*scipy.sqrt(1.-lambdaMat**2.)*jtj
    diagtermMat = (lambdaMat**2.)*jaa + (1.-lambdaMat**2.)*jbb
    return crosstermMat, diagtermMat    

def calcR2mat_3(jtj,Dmat):
#def calcR2mat_3(jtj,lambdaMat):
    n = len(jtj)
    jaa = scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj))
    r2mat = (1./2.)*(-scipy.sqrt(Dmat)+4.*(jtj**2)/scipy.sqrt(Dmat) + (jaa - 4.*jtj*scipy.sqrt((jtj**2)/Dmat) + jbb))
#    r2mat = (lambdaMat**2.)*jaa + (1.-lambdaMat**2.)*jbb + 2.*lambdaMat*scipy.sqrt(1.-lambdaMat**2.)*jtj
    return r2mat

def calcLambdaMat_3(jtj,Dmat):
    n = len(jtj)
    lambdaMat = -(1./scipy.sqrt(2.))*scipy.sqrt(1.+(-scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d')) +  scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj)))/scipy.sqrt(Dmat))
    return lambdaMat

def calcR2mat_2(jtj,Dmat):
    n = len(jtj)
    jaa = scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj))
    r2mat = (1./2.)*(scipy.sqrt(Dmat)-4.*(jtj**2)/scipy.sqrt(Dmat) + (jaa + 4.*jtj*scipy.sqrt((jtj**2)/Dmat) + jbb))
    return r2mat

def calcLambdaMat_2(jtj,Dmat):
    n = len(jtj)
    lambdaMat = (1./scipy.sqrt(2.))*scipy.sqrt(1.+(scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d')) -  scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj)))/scipy.sqrt(Dmat))
    return lambdaMat

def calcR2mat_1(jtj,Dmat):
    n = len(jtj)
    jaa = scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj))
    r2mat = (1./2.)*(scipy.sqrt(Dmat)-4.*(jtj**2)/scipy.sqrt(Dmat) + (jaa - 4.*jtj*scipy.sqrt((jtj**2)/Dmat) + jbb))
    return r2mat

def calcLambdaMat_1(jtj,Dmat):
    n = len(jtj)
    lambdaMat = -(1./scipy.sqrt(2.))*scipy.sqrt(1.+(scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d')) -  scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj)))/scipy.sqrt(Dmat))
    return lambdaMat

def calcR2mat_4(jtj,Dmat):
    n = len(jtj)
    jaa = scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d'))
    jbb = scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj))
    r2mat = (1./2.)*(-scipy.sqrt(Dmat)+4.*(jtj**2)/scipy.sqrt(Dmat) + (jaa + 4.*jtj*scipy.sqrt((jtj**2)/Dmat) + jbb))
    return r2mat

def calcLambdaMat_4(jtj,Dmat):
    n = len(jtj)
    lambdaMat = (1./scipy.sqrt(2.))*scipy.sqrt(1.+(-scipy.outer(scipy.diagonal(jtj),scipy.ones(n,'d')) +  scipy.outer(scipy.ones(n,'d'),scipy.diagonal(jtj)))/scipy.sqrt(Dmat))
    return lambdaMat

def checkLambda(jtj,lMat,r2Mat,a,b):
    lambdas = scipy.arange(-1.,1.,0.001)
    residuals = [l*l*jtj[a,a]+(1.-l*l)*jtj[b,b]+2.*l*scipy.sqrt(1.-l*l)*jtj[a,b] for l in lambdas]
    Plotting.plot(lambdas,residuals)
    Plotting.plot([lMat[a,b]],[r2Mat[a,b]],'o')
    Plotting.title(str(a)+', '+str(b))
    return

def getMinIndex(listofNums):
    minInd = 0
    minVal = listofNums[minInd]
    for ii in range(1,len(listofNums)):
        if listofNums[ii] < minVal:
            minInd = ii
            minVal = listofNums[minInd]
    return minInd

def getMinLambdaR2(lambdaMats, res2Mats):
    numMats = len(lambdaMats)
    numParams = len(lambdaMats[0])
    minR2 = scipy.zeros((numParams,numParams),'d')
    minLambda = scipy.zeros((numParams,numParams),'d')
    for ii in range(numParams):
        for jj in range(numParams):
            minIndex = getMinIndex([res2Mats[k][ii,jj] for k in range(numMats)])
            minR2[ii,jj] = res2Mats[minIndex][ii,jj]
            minLambda[ii,jj] = lambdaMats[minIndex][ii,jj]
    return minLambda, minR2

def getLambdaR2(jtj,Dmat):
    lambdaMats = [calcLambdaMat_1(jtj,Dmat),calcLambdaMat_2(jtj,Dmat), calcLambdaMat_3(jtj,Dmat),calcLambdaMat_4(jtj,Dmat)]
    r2Mats = [calcR2mat_1(jtj,Dmat),calcR2mat_2(jtj,Dmat),calcR2mat_3(jtj,Dmat),calcR2mat_4(jtj,Dmat)]
    lambdaMat, R2Mat = getMinLambdaR2(lambdaMats, r2Mats)
    return lambdaMat, R2Mat
