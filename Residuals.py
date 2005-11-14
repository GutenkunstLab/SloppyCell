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
            # Default to 0 if parameter not in senspredictions, (i.e. it may
            #  not have been in a given calculation
	    v1 = self.dy(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
	    senspredictions[self.cKey][self.yKey][self.xVal].get(pname, 0)
	    #v1 = senspredictions[self.cKey][self.yKey][self.xVal][pname]
	  else :
	    v1 = 0.0
	  if len(self.dintVar(predictions, internalVars, params)) > 0 :
	    v2 = self.dintVar(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
	    internalVarsDerivs['scaleFactors'][self.exptKey][self.yKey][pname]
	  else :
	    v2 = 0.0
	  if len(self.dp(predictions, internalVars, params)) > 0 and self.pKey == pname:
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
        return {self.pKey: 1./(params.getByKey(self.pKey)*self.sigmaLogPVal)}

    def dy(self, predictions, internalVars, params):
        return {}
    
    def dintVar(self, predictions, internalVars, params):
    	return {}

    def Dp(self,predictions,senspredictions,internalVars,internalVarsDerivs,\
    	params):
	""" Total derivative w.r.t p of the residual. Can be used as it stands 
        for other residual
	types but v3 may have to be changed depending on how it is decided 
        to key dp """
	deriv = []
	for pname in params.keys() :
            if len(self.dy(predictions, internalVars, params)) > 0 :
                # Default to 0 if parameter not in senspredictions, (i.e. it may
                #  not have been in a given calculation
                v1 = self.dy(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
                     senspredictions[self.cKey][self.yKey][self.xVal].get(pname, 0)
	        #v1 = senspredictions[self.cKey][self.yKey][self.xVal][pname]
            else :
                v1 = 0.0
            if len(self.dintVar(predictions, internalVars, params)) > 0 :
                v2 = self.dintVar(predictions, internalVars, params)[self.cKey][self.yKey][self.xVal]*\
                     internalVarsDerivs['scaleFactors'][self.exptKey][self.yKey].get(pname, 0)
            else :
                v2 = 0.0
            if len(self.dp(predictions, internalVars, params)) > 0 and self.pKey == pname:
                # dp should be keyed on the residual name? Should dy and dintVar be aswell?
                v3 = self.dp(predictions, internalVars, params).get(pname, 0)
            else :
                v3 = 0.0
            deriv.append(v1+v2+v3)
	return deriv

class PeriodCheckResidual(Residual):
    def __init__(self, key, calcKey, depVarKey, indVarValue,  depVarMeasurement,
                 depVarSigma):
        Residual.__init__(self, key)
        self.cKey = calcKey
        self.yKey = depVarKey
        self.xVal = indVarValue
        self.yMeas = depVarMeasurement
        self.ySigma = depVarSigma

    def GetRequiredVarsByCalc(self):
        # Need the points between xVal and xVal + 2*periods
        return {self.cKey: {self.yKey: [self.xVal, self.xVal+2.0*self.yMeas]}}

    def GetValue(self, predictions, internalVars, params):
        # Find the period
        traj = predictions[self.cKey][self.yKey]
        times = traj.keys()
        times.sort()

        maximums=[]
        for index in range(1,len(times)-1):
            if times[index]<self.xVal or (self.xVal+2.0*self.yMeas)<times[index]: continue
            t1,t2,t3 = times[index-1:index+2]
            y1,y2,y3 = traj[t1], traj[t2], traj[t3]
            if y1<y2 and y3<y2: # Use a quadratic approximation to find the maximum
                (a,b,c)=scipy.dot(scipy.linalg.inv([[t1**2,t1,1],[t2**2,t2,1],[t3**2,t3,1]]),[y1,y2,y3])
                maximums.append(-b/2/a)

        if len(maximums)<2: theoryVal = 2*self.yMeas
        else: theoryVal = maximums[1] - maximums[0]

        return (theoryVal - self.yMeas)/self.ySigma

    def dp(self, predictions, internalVars, params):
        return {}

    def dy(self, predictions, internalVars, params):
        return {}

    def Dp(self,predictions,senspredictions,internalVars,internalVarsDerivs,params):
	""" Total derivative w.r.t p of the residual. Can be used as it stands 
        for other residual types but v3 may have to be changed depending on how
        it is decided to key dp """
        raise NotImplementedError


class AmplitudeCheckResidual(Residual):
    def __init__(self, key, calcKey, depVarKey, indVarValue0, indVarValue1, period,
                 depVarSigma, exptKey):
        Residual.__init__(self, key)
        self.cKey = calcKey
        self.yKey = depVarKey
        self.xVal = indVarValue0
        self.xTestVal = indVarValue1
        self.period = period
        self.ySigma = depVarSigma
        self.exptKey = exptKey

    def GetRequiredVarsByCalc(self):
        # Need the points between xVal and xVal + 2*periods
        return {self.cKey: {self.yKey: [self.xVal, self.xVal+self.period, self.xTestVal, self.xTestVal+self.period]}}

    def GetValue(self, predictions, internalVars, params):
        times = predictions[self.cKey][self.yKey].keys()
        if self.exptKey in internalVars['scaleFactors'].keys() \
               and self.yKey in internalVars['scaleFactors'][self.exptKey].keys():
            scale_factor = internalVars['scaleFactors'][self.exptKey][self.yKey]
        else:
            scale_factor = 1.0

        # Get the indices of the points to use and integrate the areas
        times = predictions[self.cKey][self.yKey].keys()
        times.sort()
        startIndex,endStartIndex = times.index(self.xVal), times.index(self.xVal+self.period)
        testIndex,endTestIndex = times.index(self.xTestVal), times.index(self.xTestVal+self.period)
        
        x,y=[],[]
        for t in times[startIndex:endStartIndex+1]:
            x.append(t)
            y.append(scale_factor*predictions[self.cKey][self.yKey][t])
        measVal = scipy.integrate.simps(y,x)
            
        x,y=[],[]
        for t in times[testIndex:endTestIndex+1]:
            x.append(t)
            y.append(scale_factor*predictions[self.cKey][self.yKey][t])
        theoryVal = scipy.integrate.simps(y,x)

        return (theoryVal - measVal)/self.ySigma

    def dp(self, predictions, internalVars, params):
        return {}

    def dy(self, predictions, internalVars, params):
        return {}

    def Dp(self,predictions,senspredictions,internalVars,internalVarsDerivs,params):
	""" Total derivative w.r.t p of the residual. Can be used as it stands 
        for other residual types but v3 may have to be changed depending on how
        it is decided to key dp """
        raise NotImplementedError
