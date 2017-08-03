import scipy

class Residual:
    def __init__(self, key):
        self.key = key

    def GetKey(self):
        """
        A unique (and hashable) identifier for this residual.
        """
        return key

    def GetValue(self, predictions, internalVars, params):
        """
        The value of the residual give a current set of parameters and
        resulting predictions.
        """
        raise NotImplementedError

    def GetRequiredVarsByCalc(self):
        """
        The variables that need to be calculated to evaluate this residual.

        Should return a nested dictionary of the form:
            {calc name: {dependent var name: [independent var values]}}
        """
        raise NotImplementedError

    def dp(self, predictions, internalVars, params):
        """
        Partial derivative of the residual with respect to any parameters.

        Should return a dictionary of the form:
            {parameter name: derivative}
        """
        raise NotImplementedError

    def dy(self, predictions, internalVars, params):
        """
        Partial derivative of the residual with respect to any calculated 
        variables.

        Should return a dictionary of the form:
            {calculation name: {variable name: {x value: deriv}}}
        """
        raise NotImplementedError

    def dintVars(self, predictions, internalVars, params):
        """
        Partial derivative of the residual with respect to any internal
        variables.

        Should return a dictionary of the form:
            {type of internal var: {expt name: {variable name: derivative}}}

        XXX: This form of nesting is only appropriate for scale factors.
        """
        raise NotImplementedError

    def Dp(self, predictions, senspredictions, internalVars, internalVarsDerivs,
           params):
	""" 
        Total derivatives with respect to all parameters of the residual.

        Should return a list with the derivatives in the same order as params.

        XXX: This only works with internvalVars that are indexed like scale
        factors.
        """
	derivs_wrt_p = []
	for pname in params.keys():
            deriv = 0

            # This first term is dres/dy * dy/dp
            dres_dy = self.dy(predictions, internalVars, params)
            for calcKey in dres_dy.keys():
                for yKey in dres_dy[calcKey].keys():
                    for xVal in dres_dy[calcKey][yKey].keys():
                        dres_dy_this = dres_dy[calcKey][yKey][xVal]
                        # We default to 0 if parameter not in senspredictions. 
                        # (It may not be involved in a given calculation.)
                        dy_dp_this = senspredictions[calcKey][yKey][xVal].get(pname, 0)
                        deriv += dres_dy_this * dy_dp_this

            # This term is dres/dintVar * dintVar/dp
            dres_dintVars = self.dintVars(predictions, internalVars, params)
            for intVar_type in dres_dintVars.keys():
                for exptKey in dres_dintVars[intVar_type].keys():
                    for yKey in dres_dintVars[intVar_type][exptKey].keys():
                        dres_dintVar_this = dres_dintVars[intVar_type][exptKey][yKey]
                        dintVar_dp_this = internalVarsDerivs[intVar_type][exptKey][yKey][pname]
                        deriv += dres_dintVar_this * dintVar_dp_this

            # Finally dres/dp
            dres_dp = self.dp(predictions, internalVars, params)
            deriv += dres_dp.get(pname, 0)

            derivs_wrt_p.append(deriv)

	return derivs_wrt_p

class ScaledErrorInFit(Residual):
    def __init__(self, key, depVarKey, calcKey, indVarValue,  depVarMeasurement,
                 depVarSigma, exptKey):
        Residual.__init__(self, key)
        self.yKey = depVarKey
        self.calcKey = calcKey
        self.xVal = indVarValue
        self.yMeas = depVarMeasurement
        self.ySigma = depVarSigma
        self.exptKey = exptKey

    def GetRequiredVarsByCalc(self):
        return {self.calcKey: {self.yKey: [self.xVal]}}

    def GetValue(self, predictions, internalVars, params):
        scale_factor = internalVars['scaleFactors'][self.exptKey][self.yKey]
        raw_pred_val = predictions[self.calcKey][self.yKey][self.xVal]
        return (scale_factor * raw_pred_val - self.yMeas)/self.ySigma

    def dp(self, predictions, internalVars, params):
        return {}

    def dy(self, predictions, internalVars, params):
        scale_factor = internalVars['scaleFactors'][self.exptKey][self.yKey]
        deriv = scale_factor / self.ySigma
        return {self.calcKey: {self.yKey: {self.xVal: deriv}}}

    def dintVars(self, predictions, internalVars, params):
	raw_pred_val = predictions[self.calcKey][self.yKey][self.xVal]
        deriv = raw_pred_val / self.ySigma
        return {'scaleFactors': {self.exptKey: {self.yKey: deriv}}}

class PriorInLog(Residual):
    def __init__(self, key, pKey, logPVal, sigmaLogPVal):
        Residual.__init__(self, key)
        self.pKey = pKey
        self.logPVal = logPVal
        self.sigmaLogPVal = sigmaLogPVal

    def GetValue(self, predictions, internalVars, params):
        return (scipy.log(params.get(self.pKey)) - self.logPVal) / self.sigmaLogPVal

    def dp(self, predictions, internalVars, params):
        return {self.pKey: 1./(params.get(self.pKey) * self.sigmaLogPVal)}

    def dy(self, predictions, internalVars, params):
        return {}

    def dintVars(self, predictions, internalVars, params):
        return {}

class Prior(Residual):
    def __init__(self, key, pKey, pVal, sigmaPVal):
        Residual.__init__(self, key)
        self.pKey = pKey
        self.pVal = pVal
        self.sigmaPVal = sigmaPVal

    def GetValue(self, predictions, internalVars, params):
        return (params.get(self.pKey) - self.pVal) / self.sigmaPVal

    def dp(self, predictions, internalVars, params):
        return {self.pKey: 1./self.sigmaPVal}

    def dy(self, predictions, internalVars, params):
        return {}

    def dintVars(self, predictions, internalVars, params):
        return {}

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
            if times[index]<self.xVal\
               or (self.xVal+2.0*self.yMeas)<times[index]: continue
            t1,t2,t3 = times[index-1:index+2]
            y1,y2,y3 = traj[t1], traj[t2], traj[t3]
            if y1<y2 and y3<y2: 
                # Use a quadratic approximation to find the maximum
                (a,b,c)=scipy.dot(scipy.linalg.inv([[t1**2,t1,1],
                                                    [t2**2,t2,1],
                                                    [t3**2,t3,1]]),[y1,y2,y3])
                maximums.append(-b/2/a)

        if len(maximums)<2: 
            theoryVal = 2*self.yMeas
        else: 
            theoryVal = maximums[1] - maximums[0]

        return (theoryVal - self.yMeas)/self.ySigma

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

class IntegralDataResidual(Residual):
    def __init__(self, name, var, exptKey, calc, traj, uncert_traj, interval):
        self.name = name
        self.var = var
        self.exptKey = exptKey
        self.calc = calc
        self.traj = traj
        self.uncert_traj = uncert_traj
        self.interval = interval

    def GetValue(self, predictions, internalVars, params):
        sf = internalVars['scaleFactors'][self.exptKey][self.var]
        data_traj = self.traj
        uncert_traj = self.uncert_traj
        theory_traj = predictions[self.calc]['full trajectory']
        var = self.var
        def integrand(t):
            theory = theory_traj.evaluate_interpolated_traj(var, t)
            data = data_traj.evaluate_interpolated_traj(var, t)
            uncert = uncert_traj.evaluate_interpolated_traj(var, t)
            return (sf*theory - data)**2/uncert**2
        val, error = scipy.integrate.quad(integrand,
                                          self.interval[0], self.interval[1],
                                          limit = int(1e5))
        T = self.interval[1] - self.interval[0]
        return scipy.sqrt(val/T)

class ScaledExtremum(Residual):
    def __init__(self, key, var, calcKey, val,
                 sigma, exptKey, minTime=None, maxTime=None, type=None):
        Residual.__init__(self, key)
        self.var = var
        self.calcKey = calcKey
        self.yMeas = val
        self.ySigma = sigma
        self.exptKey = exptKey
        self.minTime,self.maxTime = minTime,maxTime
        self.last_time_result = None
        self.type = type
        if self.type == 'max':
            self.yKey= self.var + '_maximum'
        elif self.type == 'min':
            self.yKey = self.var + '_minimum'

    def GetRequiredVarsByCalc(self):
        return {self.calcKey: {self.yKey: [self.minTime, self.maxTime]}}

    def GetValue(self, predictions, internalVars, params):
        scale_factor = internalVars['scaleFactors'][self.exptKey][self.var]
        # predictions entry is time,value
        # We store the last time result for use in plotting.
        self.last_time_result, raw_pred_val = \
                predictions[self.calcKey][self.yKey][self.minTime,self.maxTime]
        return (scale_factor * raw_pred_val - self.yMeas)/self.ySigma

    def Dp(self, predictions, senspredictions, internalVars, internalVarsDerivs,
           params):
	""" 
        Total derivatives with respect to all parameters of the residual.

        Should return a list with the derivatives in the same order as params.
        """
	derivs_wrt_p = []
	for pname in params.keys():
            deriv = 0

            # This first term is dres/dy * dy/dp
            scale_factor = internalVars['scaleFactors'][self.exptKey][self.var]
            dres_dy = scale_factor / self.ySigma
            dy_dp = senspredictions[self.calcKey][self.yKey][self.minTime,self.maxTime].get(pname, 0)
            deriv += dres_dy * dy_dp

            # This term is dres/dscale_factor* dscale_factor/dp
            raw_pred_val = predictions[self.calcKey][self.yKey][self.minTime,self.maxTime][1]
            dres_dsf = raw_pred_val / self.ySigma
            dsf_dp = internalVarsDerivs['scaleFactors'][self.exptKey][self.var][pname]
            deriv += dres_dsf * dsf_dp

            derivs_wrt_p.append(deriv)

	return derivs_wrt_p
