import scipy

class Residual:
    def __init__(self, key):
        self.key = key

    def GetKey(self):
        return key

    def GetValue(self, predictions, internalVars, params):
        raise NotImplementedError

    def GetRequiredVarsByCalc(self):
        return {}

    def dp(self, predictions, internalVars, params):
        return {}

    def dy(self, predictions, internalVars, params):
        return {}

class ScaledErrorInFit(Residual):
    def __init__(self, key, depVarKey, calcKey, indVarValue,  depVarMeasurement,
                 depVarSigma, exptKey):
        Residual.__init__(self, key)
        self.yKey = depVarKey
        self.cKey = calcKey
        self.xVal = indVarValue
        self.yMeas = depVarMeasurement
        self.ySigma = depVarSigma
        self.exptKey = exptKey

    def GetRequiredVarsByCalc(self):
        return {self.cKey: {self.yKey: [self.xVal]}}

    def GetValue(self, predictions, internalVars, params):
        theoryVal = internalVars['scaleFactors'][self.exptKey][self.yKey]\
                * predictions[self.cKey][self.yKey][self.xVal]
        return (theoryVal - self.yMeas)/self.ySigma

    def dp(self, predictions, internalVars, params):
        return {}

    def dy(self, predictions, internalVars, params):
        deriv = internalVars['scaleFactors'][self.exptKey][self.yKey]\
                / self.ySigma
        return {self.cKey: {self.yKey: {self.xVal: deriv}}}
    def dintVar(self, predictions, internalVars, params):
	deriv = predictions[self.cKey][self.yKey][self.xVal]/self.ySigma
	return {self.cKey: {self.yKey: {self.xVal: deriv}}}

    def Dp(self,predictions,senspredictions,internalVars,internalVarsDerivs,\
    	params):
	""" Total derivative w.r.t p of the residual. Can be used as it stands 
        for other residual
	types but v3 may have to be changed depending on how it is decided 
        to key dp """
	deriv = []
	for pname in params.keys() :
	  if len(self.dy(predictions, internalVars, params)) > 0 :
	    v1 = self.dy(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
	    senspredictions[self.cKey][self.yKey][self.xVal][pname]
	    #v1 = senspredictions[self.cKey][self.yKey][self.xVal][pname]
	  else :
	    v1 = 0.0
	  if len(self.dintVar(predictions, internalVars, params)) > 0 :
	    v2 = self.dintVar(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
	    internalVarsDerivs['scaleFactors'][self.exptKey][self.yKey][pname]
	  else :
	    v2 = 0.0
	  if len(self.dp(predictions, internalVars, params)) > 0 and self.pKey is pname:
	    # dp should be keyed on the residual name? Should dy and dintVar be aswell?
	    v3 = self.dp(predictions, internalVars, params)[self.key]
	  else :
	    v3 = 0.0
	  deriv.append(v1+v2+v3)

	return deriv

class PriorInLog(Residual):
    def __init__(self, key, pKey, logPVal, sigmaLogPVal):
        Residual.__init__(self, key)
        self.pKey = pKey
        self.logPVal = logPVal
        self.sigmaLogPVal = sigmaLogPVal

    def GetRequiredVarsByCalc(self):
        return {}

    def GetValue(self, predictions, internalVars, params):
        return (scipy.log(params.getByKey(self.pKey)) - self.logPVal)\
                / self.sigmaLogPVal

    def dp(self, predictions, internalVars, params):
        return {pKey: scipy.log(params.getByKey(self.pKey))/self.sigmaLogPVal}

    def dy(self, predictions, internalVars, params):
        return {}
    
    def dintVar(self, predictions, internalVars, params):
    	return {}
