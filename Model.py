import copy
import exceptions
import sets
import os

import scipy

import SloppyCell
import Residuals

class Model:
    """
    Model(experiment name list, calculation name list)

    A Model unites a CalculationCollection and an ExperimentCollection. Thus
    it can calculate a cost.

    The constructor takes a list of Experiments and Calculations to include.

    XXX: In the future this class should also hold interfaces to optimization
         routines and ensemble generation. Basically anything that will work
         independent of the form of the Calculations and Experiments.
    """

    def __init__(self, exptColl, calcColl):
        self.calcVals = {}
        self.calcSensitivityVals = {}
	self.internalVars = {}
        self.internalVarsDerivs = {}
	self.residuals = SloppyCell.KeyedList()

        self.SetExperimentCollection(exptColl)
        self.SetCalculationCollection(calcColl)
	
    def CostPerResidualFromLogParams(self, params):
        return self.Cost(scipy.exp(params))/float(len(self.residuals))

    def CostFromLogParams(self, params):
    	c = self.Cost(scipy.exp(params))
	return c

    def cost(self, params):
        """
        Cost(parmeters) -> float

        Returns the sum of the squares of the residuals from
        ComputeResidualsWithScaleFactors.
        """
        self.params.update(params)
        self.CalculateForAllDataPoints(params)
        self.ComputeInternalVariables()

        resTotal = 0
        for res in self.residuals.values():
            resTotal += res.GetValue(self.calcVals, self.internalVars,
                                     self.params)**2

        return resTotal

    Cost = cost

    def AddResidual(self, res):
        self.residuals.setByKey(res.key, res)

    def Force(self, params, epsf):
        """
        Force(parameters, epsilon factor) -> list

        Returns a list containing the numerical gradient of the Cost with
        respect to each parameter (in the parameter order of the
        CalculationCollection). Each element of the gradient is:
            Cost(param + eps factor * param) - Cost(param - eps factor * param)/
              (2 * eps factor * param).
        """

        force = []
        params = scipy.array(params)
        for index, param in enumerate(params):
            eps = epsf * param

            paramsPlus = params.copy()
            paramsPlus[index] = param + eps
            costPlus = self.Cost(paramsPlus)

            paramsMinus = params.copy()
            paramsMinus[index] = param - eps
            costMinus = self.Cost(paramsMinus)

            force.append((costPlus-costMinus)/(2.0*eps))

        return force

    def CalculateForAllDataPoints(self, params):
        """
        CalculateForAllDataPoints(parameters) -> dictionary

        Gets a dictionary of measured independent variables indexed by 
        calculation from the ExperimentCollection and passes it to the
        CalculationCollection. The returned dictionary is of the form:
         dictionary[experiment][calculation][dependent variable]
                   [independent variabled] -> calculated value.
        """
        varsByCalc = self.GetExperimentCollection().GetVarsByCalc()
	self.GetCalculationCollection().Calculate(varsByCalc, params)
        self.calcVals = self.GetCalculationCollection().\
                GetResults(varsByCalc)
        return self.calcVals

    def CalculateSensitivitiesForAllDataPoints(self, params):
        """
        CalculateSensitivitiesForAllDataPoints(parameters) -> dictionary

        Gets a dictionary of measured independent variables indexed by
        calculation from the ExperimentCollection and passes it to the
        CalculationCollection. The returned dictionary is of the form:
         dictionary[experiment][calculation][dependent variable]
                   [independent variabled][parameter] -> sensitivity.
        """
        varsByCalc = self.GetExperimentCollection().GetVarsByCalc()
	self.GetCalculationCollection().CalculateSensitivity(varsByCalc, params)
        self.calcSensitivityVals = self.GetCalculationCollection().\
                GetSensitivityResults(varsByCalc)
        # might as well fill up the values for the trajectory (calcVals) 
        #  while we're at it:
	# no need to Calculate because the values are stored after the call to
        #  CalculateSensitivity
	# self.GetCalculationCollection().Calculate(varsByCalc, params)
        self.calcVals = self.GetCalculationCollection().GetResults(varsByCalc)
	return self.calcSensitivityVals

    def ComputeResidualsWithScaleFactors(self, params):
    #    """
    #    ComputeResidualsWithScaleFactors(parameters) -> dictionary

    #    Computes a dictionary of residuals corresponding to the passed in
    #    parameters. The dictionary returned is of the form:
    #      1) dictionary[(experiment, calculation, dependent variable,
    #                     independent variable)] -> residual
    #         for residuals corresponding to experimental data
    #      2) dictionary[parameter] -> residual
    #         for residuals corresponding to priors

    #    XXX: Currently only calculates residuals as
    #         (scale factor * prediction - measurement)/uncertainty
    #    XXX: The structure of the return value may need to change to accomodate
    #         more general forms of residuals.
    #    """
        self.params.update(params)
	self.CalculateSensitivitiesForAllDataPoints(params)
	#self.CalculateForAllDataPoints(params)
	self.ComputeInternalVariables()
        residuals = {}
        dataColl = self.GetExperimentCollection().GetData()
        for exptName in dataColl:
            for calcName in dataColl[exptName]:
                for depVar in dataColl[exptName][calcName]:
                    for indVar, (data, error)\
                            in dataColl[exptName][calcName][depVar].items():
                        res = (self.internalVars['scaleFactors'][exptName][depVar]\
                               *self.calcVals[calcName][depVar][indVar]-data)/error
                        residuals[(exptName, calcName, depVar, indVar)] = res

        return residuals


    def ComputeInternalVariables(self):
        self.internalVars['scaleFactors'] = self.compute_scale_factors()

    def compute_scale_factors(self):
        """
        Compute the scale factors for the current parameters and return a dict.

        The dictionary is of the form dict[exptId][varId] = scale_factor
        """

        scale_factors = {}
        for exptId, expt in self.GetExperimentCollection().items():
            scale_factors[exptId] = self._compute_sf_for_expt(expt)

        return scale_factors

    def _compute_sf_for_expt(self, expt):
        # Compute the scale factors for a given experiment
        scale_factors = {}

        exptData = expt.GetData()
        fixed_sf = expt.GetFixedScaleFactors()

        # Get the variables measured in this experiment
        measuredVars = sets.Set()
        for calcId in exptData:
            measuredVars.union_update(exptData[calcId].keys())

        sf_groups = map(sets.Set, expt.get_shared_scale_factors())
        # Flatten out the list of shared scale factors.
        flattened = []
        for g in sf_groups:
            flattened.extend(g)
        # These are variables that don't share scale factors
        unshared = [sets.Set([var]) for var in measuredVars
                    if var not in flattened]
        sf_groups.extend(unshared)

        for group in sf_groups:
            # Do any of the variables in this group have fixed scale factors?
            fixed = group.intersection(fixed_sf.keys())
            fixedAt = sets.Set([fixed_sf[var] for var in fixed])
            if len(fixedAt) == 1:
                value = fixedAt.pop()
                for var in group:
                    scale_factors[var] = value
                continue
            elif len(fixedAt) > 1:
                    raise exceptions.ValueError('Shared scale factors fixed at inconsistent values in experiment %s!' % expt.GetName())

            # Finally, compute the scale factor for this group
            theoryDotData, theoryDotTheory = 0, 0
            for calc in exptData:
                # Pull out the vars we have measured for this calculation
                for var in group.intersection(exptData[calc].keys()):
                    for indVar, (data, error) in exptData[calc][var].items():
                        theory = self.calcVals[calc][var][indVar]
                        theoryDotData += (theory * data) / error**2
                        theoryDotTheory += theory**2 / error**2

            for var in group:
                if theoryDotTheory is not 0:
                    scale_factors[var] = theoryDotData/theoryDotTheory
                else:
                    scale_factors[var] = 1

        return scale_factors

    def ComputeInternalVariableDerivs(self):
        """
        compute_scale_factorsDerivs() -> dictionary

        Returns the scale factor derivatives w.r.t. parameters
	appropriate for each chemical in each
        experiment, given the current parameters. The returned dictionary
        is of the form: internalVarsDerivs['scaleFactors'] = dict[experiment][chemical][parametername]
	 -> derivative.
        """

        self.internalVarsDerivs['scaleFactors'] = {}
	p = self.GetCalculationCollection().GetParameters()

	for exptName, expt in self.GetExperimentCollection().items():
            self.internalVarsDerivs['scaleFactors'][exptName] = {}
	    exptData = expt.GetData()

            # Get the dependent variables measured in this experiment
            exptDepVars = sets.Set()
            for calc in exptData:
                exptDepVars.union_update(expt.GetData()[calc].keys())

            for depVar in exptDepVars:
		self.internalVarsDerivs['scaleFactors'][exptName][depVar] = {}
		if depVar in expt.GetFixedScaleFactors():
                    for pname in p.keys() :
		    	self.internalVarsDerivs['scaleFactors'][exptName]\
                                [depVar][pname] = 0.0
                    continue

                theoryDotData, theoryDotTheory = 0, 0
                for calc in exptData:
                    if depVar in exptData[calc].keys():
                        for indVar, (data, error)\
                                in exptData[calc][depVar].items():
                            theory = self.calcVals[calc][depVar][indVar]
                            theoryDotData += (theory * data) / error**2
                            theoryDotTheory += theory**2 / error**2

                # now get derivative of the scalefactor
                for pname in p.keys() :
                    theoryDotData, theoryDotTheory = 0, 0

                    for calc in exptData:
                       clc = self.calcColl[calc]
                       if depVar in exptData[calc].keys():
                            for indVar, (data, error)\
                                    in exptData[calc][depVar].items():
                                theory = self.calcVals[calc][depVar][indVar]

                                theorysens = self.calcSensitivityVals[calc][depVar][indVar][pname]
                                theorysensDotData += (theorysens * data) / error**2
                                theorysensDotTheory += (theorysens * theory) / error**2

                    try:
                        self.internalVarsDerivs['scaleFactors'][exptName][depVar]\
                                [pname] = theorysensDotData/theoryDotTheory - \
                                2*theoryDotData*theorysensDotTheory/(theoryDotTheory)**2
                    except ZeroDivisionError:
                        self.internalVarsDerivs['scaleFactors'][exptName][depVar][pname] = 0

	return self.internalVarsDerivs['scaleFactors']

    def GetJacobian(self,params) :
    	"""
	GetJacobian(parameters) -> dictionary

	Gets a dictionary of the sensitivities at the time points of
	the independent variables for the measured dependent variables
	for each calculation and experiment.
	Form:
	dictionary[(experiment,calculation,dependent variable,
	independent variable)] -> result

	result is a vector of length number of parameters containing
	the sensitivity at that time point, in the order of the ordered
	parameters

	"""
	deriv = {}

	self.CalculateSensitivitiesForAllDataPoints(params)
	self.ComputeInternalVariables()
	self.ComputeInternalVariableDerivs()

	for resName, resInstance in self.residuals.items() :
	   deriv[resName] = resInstance.Dp(self.calcVals, self.calcSensitivityVals,
                                    self.internalVars, self.internalVarsDerivs, params)
	return deriv
    
    def Jacobian(self, params, epsf, relativescale=True) :
    	"""
	Finite difference the residual dictionary to get a dictionary
	for the Jacobian. It will be indexed the same as the residuals.
	Note: epsf is either a scalar or an array. If it's a scalar then
	the step sizes for each parameter are either set to this value
	or scaled by this value. If it's a vector then each parameter
	has it's own stepsize or scale for stepsize
	"""
	origParams = params.__copy__()
	j = {} # will hold the result

	res = self.ComputeResidualsWithScaleFactors(params)

	for ks in res.keys() :
		j[ks] = []

	params = scipy.array(params)

	if relativescale is True :
	       eps = epsf * abs(params)
	       for i in range(0,len(eps)) :
		  if eps[i] < min(scipy.asarray(epsf)) :
			eps[i] = min(scipy.asarray(epsf))
	else :
	       eps = epsf * scipy.ones(len(params),scipy.Float)


	for index in range(0,len(params)) :
            paramsPlus = params.copy()
            paramsPlus[index] = params[index] + eps[index]

	    resPlus = self.ComputeResidualsWithScaleFactors(paramsPlus)

	    for ks in res.keys() :
		j[ks].append((resPlus[ks]-res[ks])/(eps[index]))
	
	# NOTE: after call to ComputeResidualsWithScaleFactors the Model's
	# parameters get updated, must reset this:
	self.params.update(origParams)
	return j

    def GetJandJtJ(self,params) :

	j = self.GetJacobian(params)
	mn = scipy.zeros((len(params),len(params)),scipy.Float)

	for paramind in range(0,len(params)) :
	  for paramind1 in range(0,len(params)) :
	    sum = 0.0
	    for kys in j.keys() :
	      sum = sum + j[kys][paramind]*j[kys][paramind1]

	    mn[paramind][paramind1] = sum
    	return j,mn
   
    def GetJandJtJInLogParameters(self,params) :
    	# Formula below is exact if you have perfect data. If you don't
	# have perfect data (residuals != 0) you get an extra term when you
	# compute d^2(Cost)/(dlogp[i]dlogp[j]) which is
	# sum_resname (residual[resname] * jac[resname][j] * delta_jk * p[k])
	# but can be ignored when residuals are zeros, and maybe should be
	# ignored altogether because it can make the Hessian non-positive
	# definite
	jac, jtj = self.GetJandJtJ(params)
	for i in range(len(params)) :
            for j in range(len(params)) :
                jtj[i][j] = jtj[i][j]*params[i]*params[j]
	# extra term
	residuals = self.ComputeResidualsWithScaleFactors(params)
	for resname in residuals.keys() :
            for j in range(len(params)) :
                jtj[j][j] += residuals[resname]*jac[resname][j]*params[j]
                jac[resname][j] = jac[resname][j]*params[j]

	return jac,jtj

    def ComputeHessianElement(self, costFunc, chiSq, 
                              params, i, j, epsi, epsj, 
                              relativeScale, stepSizeCutoff):
        origPi, origPj = params[i], params[j]

        if relativeScale:
            # Steps sizes are given by eps*the value of the parameter,
            #  but the minimum step size is stepSizeCutoff
            hi, hj = scipy.maximum((epsi*abs(origPi), epsj*abs(origPj)), 
                                   (stepSizeCutoff, stepSizeCutoff))
        else:
            hi, hj = epsi, epsj

        if i == j:
            params[i] = origPi + hi
            chiSqp = costFunc(params)
            params[i] = origPi - hi
            chiSqm = costFunc(params)

            element = 0.5*(chiSqp - 2*chiSq + chiSqm)/hi**2
        else:
            ## f(xi + hi, xj + h)
            params[i] = origPi + hi
            params[j] = origPj + hj
            chiSqpp = costFunc(params)

            ## f(xi + hi, xj - hj)
            params[i] = origPi + hi
            params[j] = origPj - hj
            chiSqpm = costFunc(params)

            ## f(xi - hi, xj + hj)
            params[i] = origPi - hi
            params[j] = origPj + hj
            chiSqmp = costFunc(params)

            ## f(xi - hi, xj - hj)
            params[i] = origPi - hi
            params[j] = origPj - hj
            chiSqmm = costFunc(params)

            element = 0.5*(chiSqpp - chiSqpm - chiSqmp + chiSqmm)/(4 * hi * hj)

        params[i], params[j] = origPi, origPj

        return element


    def CalcHessianInLogParameters(self, params, eps, relativeScale = False, 
                                   stepSizeCutoff = 1e-6):
	nOv = len(params)
        if len(scipy.asarray(eps)) == 1:
            eps = scipy.ones(len(params), scipy.Float) * eps

	## compute residuals/cost at f(x)
	chiSq = self.CostFromLogParams(scipy.log(params))

	hess = scipy.zeros((nOv, nOv), scipy.Float)

	## compute all (numParams*(numParams + 1))/2 unique hessian elements
        for i in range(nOv):
            for j in range(i, nOv):
                hess[i][j] = self.ComputeHessianElement(self.CostFromLogParams,
                                                        chiSq,
                                                        scipy.log(params), 
                                                        i, j, eps[i], eps[j], 
                                                        relativeScale, 
                                                        stepSizeCutoff)
                hess[j][i] = hess[i][j]

        return hess

    def CalcHessian(self, params, eps, relativeScale = True, 
                    stepSizeCutoff = 1e-6):
	nOv = len(params)
        if len(scipy.asarray(eps)) == 1:
            eps = scipy.ones(len(params), scipy.Float) * eps

	## compute residuals/cost at f(x)
	chiSq = self.Cost(params)

	hess = scipy.zeros((nOv, nOv), scipy.Float)

	## compute all (numParams*(numParams + 1))/2 unique hessian elements
        for i in range(nOv):
            for j in range(i, nOv):
                hess[i][j] = self.ComputeHessianElement(self.Cost, chiSq,
                                                        params, i, j, 
                                                        eps[i], eps[j], 
                                                        relativeScale, 
                                                        stepSizeCutoff)
                hess[j][i] = hess[i][j]

        return hess

    def CalcHessianUsingResiduals(self,params,epsf,relativescale = True) :
    	currParams = copy.copy(params)
	nOp = len(currParams)
	J,jtj = self.GetJandJtJ(currParams)
	self.CalculateForAllDataPoints(currParams)
	self.ComputeInternalVariables()
	localCalcVals = self.calcVals
	localIntVars = self.internalVars

	paramlist = scipy.array(currParams)
	if relativescale is True :
	       eps = epsf * abs(paramlist)
	       for i in range(0,len(eps)) :
		  if eps[i] < min(scipy.asarray(epsf)) :
			eps[i] = min(scipy.asarray(epsf))
	else :
	       eps = epsf * scipy.ones(len(paramlist),scipy.Float)

	secondDeriv = scipy.zeros((nOp,nOp),scipy.Float)
	for index in range(0,len(params)) :
            paramsPlus = currParams.__copy__()
	    paramsPlus[index] = paramsPlus[index] + eps[index]
	    JPlus = self.GetJacobian(paramsPlus)

	    for res in J.keys() :
		resvalue = self.residuals.getByKey(res).GetValue\
					(localCalcVals, localIntVars, currParams)
		secondDeriv[index,:] += (scipy.array(JPlus[res])-scipy.array(J[res]))/ \
					eps[index]*resvalue

	self.params.update(currParams)
	# force symmetry, because it will never be exactly symmetric.
	# Should check that (jtj+secondDeriv) is approximately symmetric if concerned
	# about accuracy of the calculation --- the difference between the i,j and the
	# j,i entries should give you an indication of the accuracy
	hess = 0.5*(scipy.transpose(jtj + secondDeriv) + (jtj + secondDeriv))
	return hess



    ############################################################################
    # Getting/Setting variables below

    def SetExperimentCollection(self, exptColl):
        self.exptColl = exptColl

        for exptKey, expt in exptColl.items():
            exptData = expt.GetData()
            for calcKey, calcData in exptData.items():
                for depVarKey, depVarData in calcData.items():
                    for indVar, (value, uncert) in depVarData.items():
                        resName = (exptKey, calcKey, depVarKey, indVar)
                        res = Residuals.ScaledErrorInFit(resName, depVarKey,
                                                         calcKey, indVar, value,
                                                         uncert, exptKey)
                        self.residuals.setByKey(resName, res)


    def get_expts(self):
        return self.exptColl

    GetExperimentCollection = get_expts

    def SetCalculationCollection(self, calcColl):
        self.calcColl = calcColl
        self.params = calcColl.GetParameters()

    def get_calcs(self):
        return self.calcColl

    GetCalculationCollection = get_calcs

    def GetScaleFactors(self):
        return self.internalVars['scaleFactors']

    def GetResiduals(self):
        return self.residuals

    def GetCalculatedValues(self):
        return self.calcVals

    def GetInternalVariables(self):
        return self.internalVars

