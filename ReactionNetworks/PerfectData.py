import sys

import scipy

sigmaFunc = lambda traj, dataId: 1

def calculatePerfectDataLMHessian(traj, dataIds = None, optIds = None, 
                                  fixedScaleFactors = False, 
                                  returnDict = False):
    if dataIds is None:
        dataIds = traj.dynamicVarKeys + traj.assignedVarKeys
    if optIds is None:
        optIds = traj.optimizableVarKeys

    scaleFactorDerivs, LMHessianDict = {}, {}
    for dataIndex, dataId in enumerate(dataIds):
        dataSigma = sigmaFunc(traj, dataId)
        if scipy.isscalar(dataSigma):
            dataSigma = 0*traj.getVariableTrajectory(dataId) + dataSigma

        if fixedScaleFactors:
            scaleFactorDerivs[dataId] = dict(zip(optIds, 
                                                 scipy.zeros(len(optIds))))
        else:
            scaleFactorDerivs[dataId] = \
                    calculateScaleFactorDerivs(traj, dataId, dataSigma, optIds)

        LMHessianDict[dataId] = \
                computeLMHessianContribution(traj, dataId, dataSigma, optIds,
                                             scaleFactorDerivs[dataId])

    LMHessian = scipy.sum(LMHessianDict.values())
    returns = [LMHessian]
    if returnDict:
        returns.append(LMHessianDict)

    return returns

def computeIntervals(traj):
    # We want to break up our integrals when events fire, so first we figure out
    #  when they fired by looking for duplicated times in the trajectory
    eventIndices = scipy.compress(scipy.diff(traj.timepoints) == 0, 
                                  scipy.arange(len(traj.timepoints)))
    intervals = zip([0] + list(eventIndices + 1), 
                    list(eventIndices + 1) + [len(traj.timepoints)])

    return intervals


def interpolateInterval(traj, dataId, dataSigma, startIndex, endIndex):
    k = min(5, endIndex-startIndex - 1)

    times = traj.timepoints[startIndex:endIndex]
    y = traj.getVariableTrajectory(dataId)[startIndex:endIndex]
    sigma = dataSigma[startIndex:endIndex]

    yTCK = scipy.interpolate.splrep(times, y, k = k, s = 0)
    sigmaTCK = scipy.interpolate.splrep(times, sigma, k = k, s = 0)

    return yTCK, sigmaTCK

def calculateScaleFactorDerivs(traj, dataId, dataSigma, optIds):
    scaleFactorDerivs = {}
    intTheorySq = 0
    for startIndex, endIndex in computeIntervals(traj):
        k = min(5, endIndex-startIndex - 1)
        yTCK, sigmaTCK = interpolateInterval(traj, dataId, dataSigma,
                                             startIndex, endIndex)
    
        value, error = scipy.integrate.quad(_theorySqIntegrand, 
                                            traj.timepoints[startIndex], 
                                            traj.timepoints[endIndex - 1], 
                                            args = (yTCK, sigmaTCK),
                                            limit = endIndex-startIndex)
        intTheorySq += value
    
        for optId in optIds:
            sens = traj.getVariableTrajectory((dataId, optId))\
                    [startIndex:endIndex]
            optValue = traj.constantVarValues.getByKey(optId)
            sensTCK = scipy.interpolate.splrep(traj.timepoints\
                                               [startIndex:endIndex], 
                                               sens * optValue,
                                               k = k, s = 0)
    
            numerator, error = scipy.\
                    integrate.quad(_scaleFactorDerivsIntegrand,
                                   traj.timepoints[startIndex], 
                                   traj.timepoints[endIndex - 1], 
                                   args = (sensTCK, yTCK, sigmaTCK),
                                   limit = endIndex-startIndex)
    
            scaleFactorDerivs.setdefault(optId, 0)
            scaleFactorDerivs[optId] += -numerator

    for optId in optIds:
        scaleFactorDerivs[optId] /= intTheorySq

    return scaleFactorDerivs

def computeLMHessianContribution(traj, dataId, dataSigma, optIds, 
                                 scaleFactorDerivs):
    LMHessian = scipy.zeros((len(optIds), len(optIds)), scipy.Float)
    for startIndex, endIndex in computeIntervals(traj):
        k = min(5, endIndex-startIndex - 1)
        yTCK, sigmaTCK = interpolateInterval(traj, dataId, dataSigma, 
                                             startIndex, endIndex)

        sensTCKs = {}
        for optIndex, optId in enumerate(optIds):
            sens = traj.getVariableTrajectory((dataId, optId))\
                    [startIndex:endIndex]
            optValue = traj.constantVarValues.getByKey(optId)
            sensTCK = scipy.interpolate.splrep(traj.timepoints\
                                               [startIndex:endIndex], 
                                               sens * optValue,
                                               k = k, s = 0)
            sensTCKs[optId] = sensTCK

        for optIndex1, optId1 in enumerate(optIds):
            sens1TCK = sensTCKs[optId1]
            for jj, optId2 in enumerate(optIds[optIndex1:]):
                optIndex2 = jj + optIndex1

                sens2TCK = sensTCKs[optId2]

                if scaleFactorDerivs[optId1] == 0\
                   and scaleFactorDerivs[optId2] == 0:
                    integrand = _integrandFixedSF
                else:
                    integrand = _integrandVarSF

                value, error = scipy.\
                        integrate.quad(integrand, 
                                       traj.timepoints[startIndex], 
                                       traj.timepoints[endIndex - 1], 
                                       args = (sens1TCK, sens2TCK,
                                               yTCK, sigmaTCK,
                                               scaleFactorDerivs[optId1],
                                               scaleFactorDerivs[optId2]),
                                       limit = endIndex-startIndex)

                LMHessian[optIndex1][optIndex2] += value
                if optIndex1 != optIndex2:
                    LMHessian[optIndex2][optIndex1] += value
    
    LMHessian /= (traj.timepoints[-1] - traj.timepoints[0])

    return LMHessian

def _scaleFactorDerivsIntegrand(t, sensTCK, yTCK, sigmaTCK):
    sensVal = scipy.interpolate.splev(t, sensTCK)
    y = scipy.interpolate.splev(t, yTCK)
    sigma = scipy.interpolate.splev(t, sigmaTCK)
    return sensVal * y/sigma**2

def _theorySqIntegrand(t, yTCK, sigmaTCK):
    y = scipy.interpolate.splev(t, yTCK)
    sigma = scipy.interpolate.splev(t, sigmaTCK)
    return (y/sigma)**2

def _integrandFixedSF(t, sens1TCK, sens2TCK, yTCK, sigmaTCK, dB1, dB2):
    sens1Val = scipy.interpolate.splev(t, sens1TCK)
    sens2Val = scipy.interpolate.splev(t, sens2TCK)
    sigma = scipy.interpolate.splev(t, sigmaTCK)
    return sens1Val * sens2Val / sigma**2

def _integrandVarSF(t, sens1TCK, sens2TCK, yTCK, sigmaTCK, dB1, dB2):
    sens1Val = scipy.interpolate.splev(t, sens1TCK)
    sens2Val = scipy.interpolate.splev(t, sens2TCK)
    y = scipy.interpolate.splev(t, yTCK)
    sigma = scipy.interpolate.splev(t, sigmaTCK)
    return (sens1Val + dB1 * y) * (sens2Val + dB2 * y)/ sigma**2

try:
    import psyco
    psyco.bind(_theorySqIntegrand)
    psyco.bind(_scaleFactorDerivsIntegrand)
    psyco.bind(_integrandFixedSF)
    psyco.bind(_integrandVarSF)
except ImportError:
    pass
