"""
Model class that unites theory with data.
"""

__docformat__ = "restructuredtext en"

import logging
logger = logging.getLogger('Model_mod')

import copy
import sets

import scipy
try:
    import scipy.misc.limits as limits
except ImportError:
    limits = scipy.limits

import SloppyCell.Residuals as Residuals
import SloppyCell.Collections as Collections
import KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

import os
logger.debug('cwd before importing pypar: %s' % os.getcwd())
try:
    import pypar
    HAVE_PYPAR=True
except:
    HAVE_PYPAR=False
logger.debug('cwd after importing pypar: %s' % os.getcwd())
    
class Model:
    """
    A Model object connects a set of experimental data with the objects used to
    model that data.

    Most importantly, a Model can calculate a cost for a given set of 
    parameters, characterizing how well those parameters fit the data contained
    within the model.
    """
    def __init__(self, expts, calcs):
        """
        expts  A sequence of Experiments to be fit to.
        calcs  A sequence of calculation objects referred to by the 
               Experiments.
        """

        self.calcVals = {}
        self.calcSensitivityVals = {}
	self.internalVars = {}
        self.internalVarsDerivs = {}
	self.residuals = KeyedList()

        if isinstance(expts, list):
            expts = Collections.ExperimentCollection(expts)
        elif isinstance(expts, dict):
            expts = Collections.ExperimentCollection(expts.values())
        self.SetExperimentCollection(expts)

        if isinstance(calcs, list):
            calcs = Collections.CalculationCollection(calcs)
        elif isinstance(expts, dict):
            calcs = Collections.CalculationCollection(calcs.values())
        self.SetCalculationCollection(calcs)

        # Echo the cost or chi-squared on each evaluation
        self.costVerbose = False
        self.chisqVerbose = False

        self.observers = KeyedList()

    def MasterSwitch(self, flag=-1):
        """
        This is the master switch which the master calls
        to get the workers into the function of interest.
        The flag indicates which function the workers should
        call. (Don't use negative flags, they may be truncated.)

        In a parallel-ly executed script, the format should be

            if pypar.rank()==0:
                # Doing some calculations with model 'm'.
                # Anything done outside of this loop WILL BE
                # DONE BY ALL PROCESSES.
            m.MasterSwitch()

        Thus, the workers will hit this function first with
        a flag of -1, telling them to start looping over
        model 'm' waiting for the root to give them work
        (within the if statement). When the root is done
        with the work inside the loop, it will hit this fx
        with a flag of -1, telling it to signal the workers
        to stop and go back to the original script.

        Please note that, in parallel, any call of Model functions:
            CaclulateSensitivitiesForAllDataPoints
        outside of code with this structure will likely result
        in an error, since it will confuse the communication.
        """
        if not HAVE_PYPAR: return
        
        import pypar
        myid=pypar.rank()

        # Initiate looping the workers over this model
        if flag==-1:
            if myid==0: # root node will stop the workers
                flag=0
            else:
                while True:
                    if self.MasterSwitch(0): return
                
        # Broadcast the flag from the root
        flag=pypar.broadcast(flag,0)

        if flag==0: # Exit
            return True
        elif flag==1:
            # Synchronize the parameters and initial conditions
            self.synchronize()
            
            if myid!=0: # Send the workers off to the function
                temp=self.CalculateSensitivitiesForAllDataPoints(self.params)

        return False

    def synchronize(self):
        """
        This function takes care of transfering
        over the parameters of this class (self.params)
        to all the worker nodes. 

        Also, this function synchronizes the initial
        conditions, incase the user has modified them.

        Please make sure the parameters and initial
        conditions are set in the model (i.e. self.params, etc.)
        before calling this function.
        """
        # Pickle and broadcast the parameters
        self.params = self.pickle_broadcast(self.params,0)
        
        # Pickle and broadcast the initial conditions
        self.set_ICs(self.pickle_broadcast(self.get_ICs(),0))

    def pickle_broadcast(self, data, root=0):
        """
        This function will pickle and broadcast the data,
        returning the data received.

        ** This function could be put elsewhere,
        ** it has no model dependence

        Inputs:
          data -- data to send
          root -- root node
        """
        if not HAVE_PYPAR: return

        import pypar, pickle, string
        myid = pypar.rank()

        s = pickle.dumps(data) # Pickle the data

        # Synchronize the string lengths
        l = pypar.broadcast(len(s),root)
        if myid!=root: s = string.zfill('',l)

        pypar.broadcast_string(s,root)

        return pickle.loads(s)
        
    def copy(self):
        return copy.deepcopy(self)
        
    def get_params(self):
        """
        Return a copy of the current model parameters
        """
        return self.calcColl.GetParameters()

    def get_ICs(self):
        """
        Get the initial conditions currently present in a model
        for dynamic variables that are not assigned variables.

        Outputs:
          KeyedList with keys (calcName,varName) --> initialValue
        """
        ics=KeyedList()
        for calcName, calc in self.calcColl.items():
            for varName in calc.dynamicVars.keys():
                if varName in calc.assignedVars.keys(): continue
                ics.set( (calcName,varName), calc.get_var_ic(varName))
        return ics

    def set_ICs(self, ics):
        """
        Sets the initial conditions into the model. Uses the input
        format defined by 'getICs'.

        Inputs:
         ics -- Initial conditions to set in KeyedList form:
                  keys: (calcName, varName) --> intialValue

        Outputs:
         None
        """
        for (calcName, varName), initialValue in ics.items():
            self.calcColl[calcName].set_var_ic(varName, initialValue)

    def _evaluate(self, params):
        """
        Evaluate the cost for the model, returning the intermediate residuals,
        and chi-squared.

        (Summing up the residuals is a negligible amount of work. This 
         arrangment makes notification of observers much simpler.)
        """
        self.params.update(params)
        self.CalculateForAllDataPoints(params)
        self.ComputeInternalVariables()

        resvals = [res.GetValue(self.calcVals, self.internalVars, self.params)
                   for res in self.residuals.values()]

        chisq = scipy.sum(scipy.asarray(resvals)**2)
        if scipy.isnan(chisq):
            print 'Warning: Chi-Squared in NAN, setting to INF'
            chisq = scipy.inf
        cost = 0.5 * chisq
        self._notify(event = 'evaluation', 
                     resvals = resvals,
                     chisq = chisq,
                     cost = cost, 
                     params = self.params)

        return resvals, chisq, cost


    def res(self, params):
        """
        Return the residual values of the model fit given a set of parameters
        """

        return self._evaluate(params)[0]

    def res_log_params(self, logparams):
        """
        Return the residual values given the logarithm of the parameters
        """
        return self.res(scipy.exp(logparams))

    def res_dict(self, params):
        """
        Return the residual values of the model fit given a set of parameters
        in dictionary form.
        """
        return dict(zip(self.residuals.keys(), self.res(params)))

    def chisq(self, params):
        """
        Return the sum of the squares of the residuals for the model
        """
        return self._evaluate(params)[1]

    def redchisq(self, params):
        """
        Return chi-squared divided by the number of degrees of freedom

        Question: Are priors to be included in the N data points?
                  How do scale factors change the number of d.o.f.?
        """
        return self.chisq(params)/(len(self.residuals) - len(self.params))

    def cost(self, params):
        """
        Return the cost (1/2 chisq) of the model
        """
        return self._evaluate(params)[2]

        return cost

    def cost_log_params(self, log_params):
        """
        Return the cost given the logarithm of the input parameters
        """
        return self.cost(scipy.exp(log_params))


    def _notify(self, **args):
        """
        Call all observers with the given arguments.
        """
        for obs in self.observers:
            obs(**args)

    def attach_observer(self, obs_key, observer):
        """
        Add an observer to be notified by this Model.
        """
        self.observers.set(obs_key, observer)

    def detach_observer(self, obs_key):
        """
        Remove an observer from the Model.
        """
        self.observers.remove_by_key(obs_key)

    def get_observers(self):
        """
        Return the KeyedList of observers for this model.
        """
        return self.observers

    def reset_observers(self):
        """
        Call reset() for all attached observers.
        """
        for obs in self.observers:
            if hasattr(obs, 'reset'):
                obs.reset()

	
    # Deprecating...
    def CostPerResidualFromLogParams(self, params):
        return self.Cost(scipy.exp(params))/float(len(self.residuals))

    def Cost(self, params):
        return 2*self.cost(params)

    def CostFromLogParams(self, params):
    	return self.Cost(scipy.exp(params))

    resDict = res_dict
    # ...


    def AddResidual(self, res):
        self.residuals.setByKey(res.key, res)

    def Force(self, params, epsf, relativeScale=False, stepSizeCutoff=None):
        """
        Force(parameters, epsilon factor) -> list

        Returns a list containing the numerical gradient of the Cost with
        respect to each parameter (in the parameter order of the
        CalculationCollection). Each element of the gradient is:
            Cost(param + eps) - Cost(param - eps)/(2 * eps).
        If relativeScale is False then epsf is the stepsize used (it should
        already be multiplied by typicalValues before Jacobian is called)
        If relativeScale is True then epsf is multiplied by params.
        The two previous statements hold for both scalar and vector valued
        epsf.
        """

        force = []
        params = scipy.array(params)
        
        if stepSizeCutoff==None:
            stepSizeCutoff = scipy.sqrt(limits.double_epsilon)
            
	if relativeScale is True :
            eps = epsf * abs(params)
	else :
            eps = epsf * scipy.ones(len(params),scipy.Float)

        for i in range(0,len(eps)):
            if eps[i] < stepSizeCutoff:
                eps[i] = stepSizeCutoff

        for index, param in enumerate(params):
            paramsPlus = params.copy()
            paramsPlus[index] = param + eps[index]
            costPlus = self.Cost(paramsPlus)

            paramsMinus = params.copy()
            paramsMinus[index] = param - eps[index]
            costMinus = self.Cost(paramsMinus)

            force.append((costPlus-costMinus)/(2.0*eps[index]))

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
        if HAVE_PYPAR:
            import pypar
            if pypar.rank()==0 and pypar.size()>1:
                self.params.update(params)
                self.MasterSwitch(1)

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
        fixed_sf = expt.get_fixed_sf()

        # Get the variables measured in this experiment
        measuredVars = sets.Set()
        for calcId in exptData:
            measuredVars.union_update(sets.Set(exptData[calcId].keys()))

        sf_groups = map(sets.Set, expt.get_shared_sf())
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
            fixed = group.intersection(sets.Set(fixed_sf.keys()))
            fixedAt = sets.Set([fixed_sf[var] for var in fixed])
            if len(fixedAt) == 1:
                value = fixedAt.pop()
                for var in group:
                    scale_factors[var] = value
                continue
            elif len(fixedAt) > 1:
                    raise ValueError('Shared scale factors fixed at inconsistent values in experiment %s!' % expt.GetName())

            # Finally, compute the scale factor for this group
            theoryDotData, theoryDotTheory = 0, 0
            for calc in exptData:
                # Pull out the vars we have measured for this calculation
                for var in group.intersection(sets.Set(exptData[calc].keys())):
                    for indVar, (data, error) in exptData[calc][var].items():
                        theory = self.calcVals[calc][var][indVar]
                        theoryDotData += (theory * data) / error**2
                        theoryDotTheory += theory**2 / error**2

            for var in group:
                if theoryDotTheory != 0:
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
                exptDepVars.union_update(sets.Set(expt.GetData()[calc].keys()))

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
                    theorysensDotData, theorysensDotTheory = 0, 0

                    for calc in exptData:
                       clc = self.calcColl[calc]
                       if depVar in exptData[calc].keys():
                            for indVar, (data, error)\
                                    in exptData[calc][depVar].items():
                                theory = self.calcVals[calc][depVar][indVar]

                                # Default to 0 if sensitivity not calculated for
                                #  that parameter (i.e. it's not in the 
                                #  Calculation)
                                theorysens = self.calcSensitivityVals[calc][depVar][indVar].get(pname, 0.0)
                                theorysensDotData += (theorysens * data) / error**2
                                theorysensDotTheory += (theorysens * theory) / error**2

                    try:
                        self.internalVarsDerivs['scaleFactors'][exptName][depVar]\
                                [pname] = theorysensDotData/theoryDotTheory - \
                                2*theoryDotData*theorysensDotTheory/(theoryDotTheory)**2
                    except ZeroDivisionError:
                        self.internalVarsDerivs['scaleFactors'][exptName][depVar][pname] = 0

	return self.internalVarsDerivs['scaleFactors']

    def jacobian_sens(self, params):
    	"""
	Return a dictionary of the derivatives of the model residuals w.r.t.
        parameters.

        The method uses the sensitivity integration. As such, it will only
        work with ReactionNetworks.

        The dictionary is of the form:
            dict[resId] = [dres/dp1, dres/dp2...]
	"""
        self.params.update(params)

        # Calculate sensitivities
	self.CalculateSensitivitiesForAllDataPoints(params)
	self.ComputeInternalVariables()
	self.ComputeInternalVariableDerivs()

        # Calculate residual derivatives
        deriv = dict([(resId, res.Dp(self.calcVals, self.calcSensitivityVals,
                                     self.internalVars, self.internalVarsDerivs,
                                     self.params))
                      for (resId, res) in self.residuals.items()])

	return deriv

    def jacobian_fd(self, params, eps, 
                    relativeScale=False, stepSizeCutoff=None):
    	"""
	Return a dictionary of the derivatives of the model residuals w.r.t.
        parameters.

        The method uses finite differences.

        Inputs:
         params -- Parameters about which to calculate the jacobian
         eps -- Step size to take, may be vector or scalar.
         relativeScale -- If true, the eps is taken to be the fractional
                          change in parameter to use in finite differences.
         stepSizeCutoff -- Minimum step size to take.
        """
        res = self.resDict(params)

	orig_vals = scipy.array(params)

        if stepSizeCutoff is None:
            stepSizeCutoff = scipy.sqrt(limits.double_epsilon)
            
	if relativeScale:
            eps_l = scipy.maximum(eps * abs(params), stepSizeCutoff)
	else:
            eps_l = scipy.maximum(eps * scipy.ones(len(params),scipy.Float),
                                  stepSizeCutoff)

	J = {} # will hold the result
	for resId in res.keys():
            J[resId] = []
        # Two-sided finite difference
	for ii in range(len(params)):
            params[ii] = orig_vals[ii] + eps_l[ii]
	    resPlus = self.resDict(params)

            params[ii] = orig_vals[ii] - eps_l[ii]
            resMinus = self.resDict(params)

            params[ii] = orig_vals[ii]

	    for resId in res.keys() :
                J[resId].append((resPlus[resId]-resMinus[resId])/(2.*eps_l[ii]))
	
	# NOTE: after call to ComputeResidualsWithScaleFactors the Model's
	# parameters get updated, must reset this:
        self.params.update(params)
	return J

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
        return self.jacobian_sens(params)
    
    def Jacobian(self, params, epsf, relativeScale=False, stepSizeCutoff=None):
    	"""
	Finite difference the residual dictionary to get a dictionary
	for the Jacobian. It will be indexed the same as the residuals.
	Note: epsf is either a scalar or an array.
        If relativeScale is False then epsf is the stepsize used (it should
        already be multiplied by typicalValues before Jacobian is called)
        If relativeScale is True then epsf is multiplied by params.
        The two previous statements hold for both scalar and vector valued
        epsf.
        """
        return self.jacobian_fd(params, epsf, 
                                relativeScale, stepSizeCutoff)

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
	# ignored altogether because it can make the Hessian approximation non-positive
	# definite
	pnolog = scipy.exp(params)	
	jac, jtj = self.GetJandJtJ(pnolog)
	for i in range(len(params)) :
            for j in range(len(params)) :
                jtj[i][j] = jtj[i][j]*pnolog[i]*pnolog[j]

	res = self.resDict(pnolog)	
	for resname in self.residuals.keys() :
            for j in range(len(params)) :
               	# extra term --- not including it 
		# jtj[j][j] += res[resname]*jac[resname][j]*pnolog[j]
                jac[resname][j] = jac[resname][j]*pnolog[j]

	return jac,jtj

    def hessian_elem(self, func, f0, params, i, j, epsi, epsj,
                     relativeScale, stepSizeCutoff, verbose):
        """
        Return the second partial derivative for func w.r.t. parameters i and j

        f0: The value of the function at params
        eps: Sets the stepsize to try
        relativeScale: If True, step i is of size p[i] * eps, otherwise it is
                       eps
        stepSizeCutoff: The minimum stepsize to take
        """
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
            fp = func(params)

            params[i] = origPi - hi
            fm = func(params)

            element = (fp - 2*f0 + fm)/hi**2
        else:
            ## f(xi + hi, xj + h)
            params[i] = origPi + hi
            params[j] = origPj + hj
            fpp = func(params)

            ## f(xi + hi, xj - hj)
            params[i] = origPi + hi
            params[j] = origPj - hj
            fpm = func(params)

            ## f(xi - hi, xj + hj)
            params[i] = origPi - hi
            params[j] = origPj + hj
            fmp = func(params)

            ## f(xi - hi, xj - hj)
            params[i] = origPi - hi
            params[j] = origPj - hj
            fmm = func(params)

            element = (fpp - fpm - fmp + fmm)/(4 * hi * hj)

        params[i], params[j] = origPi, origPj

        self._notify(event = 'hessian element', i = i, j = j, 
                     element = element)
        if verbose: 
            print 'hessian[%i, %i] = %g' % (i, j, element)

        return element

    def hessian(self, params, epsf, relativeScale = True, 
                stepSizeCutoff = None, jacobian = None, 
                verbose = False):
    	"""
        Returns the hessian of the model.

        epsf: Sets the stepsize to try
        relativeScale: If True, step i is of size p[i] * eps, otherwise it is
                       eps
        stepSizeCutoff: The minimum stepsize to take
        jacobian: If the jacobian is passed, it will be used to estimate
                  the step size to take.
        vebose: If True, a message will be printed with each hessian element
                calculated
        """

	nOv = len(params)
        if stepSizeCutoff is None:
            stepSizeCutoff = scipy.sqrt(limits.double_epsilon)
            
        params = scipy.asarray(params)
	if relativeScale:
            eps = epsf * abs(params)
	else:
            eps = epsf * scipy.ones(len(params),scipy.Float)

        # Make sure we don't take steps smaller than stepSizeCutoff
        eps = scipy.maximum(eps, stepSizeCutoff)

        if jacobian is not None:
            # Turn off the relative scaling since that would overwrite all this
            relativeScale = False

            jacobian = scipy.asarray(jacobian)
            if len(jacobian.shape) == 0:
                resDict = self.resDict(params)
                new_jacobian = scipy.zeros(len(params),scipy.Float)
                for key, value in resDict.items():
                    new_jacobian += 2.0*value*scipy.array(jacobian[0][key])
                jacobian = new_jacobian
            elif len(jacobian.shape) == 2: # Need to sum up the total jacobian
                residuals = scipy.asarray(self.res(params))
                jacobian = 2.0*residuals*jacobian

            # If parameters are independent, then
            #  epsilon should be (sqrt(2)*J[i])^-1
            factor = 1.0/scipy.sqrt(2)
            for i in range(nOv):
                if jacobian[i] == 0.0:
                    eps[i] = 0.5*abs(params[i])
                else: 
                    # larger than stepSizeCutoff, but not more than 
                    #  half of the original parameter value
                    eps[i] = min(max(factor/abs(jacobian[i]), stepSizeCutoff), 
                                 0.5*abs(params[i]))

	## compute cost at f(x)
	f0 = self.cost(params)

	hess = scipy.zeros((nOv, nOv), scipy.Float)

	## compute all (numParams*(numParams + 1))/2 unique hessian elements
        for i in range(nOv):
            for j in range(i, nOv):
                hess[i][j] = self.hessian_elem(self.cost, f0,
                                               params, i, j, 
                                               eps[i], eps[j], 
                                               relativeScale, stepSizeCutoff,
                                               verbose)
                hess[j][i] = hess[i][j]

        return hess

    def hessian_log_params(self, params, eps,
                           relativeScale=False, stepSizeCutoff=1e-6,
                           verbose=False):
        """
        Returns the hessian of the model in log parameters.

        eps: Sets the stepsize to try
        relativeScale: If True, step i is of size p[i] * eps, otherwise it is
                       eps
        stepSizeCutoff: The minimum stepsize to take
        vebose: If True, a message will be printed with each hessian element
                calculated
        """
	nOv = len(params)
        if scipy.isscalar(eps):
            eps = scipy.ones(len(params), scipy.Float) * eps

	## compute cost at f(x)
	f0 = self.cost_log_params(scipy.log(params))

	hess = scipy.zeros((nOv, nOv), scipy.Float)

	## compute all (numParams*(numParams + 1))/2 unique hessian elements
        for i in range(nOv):
            for j in range(i, nOv):
                hess[i][j] = self.hessian_elem(self.cost_log_params, f0,
                                               scipy.log(params), 
                                               i, j, eps[i], eps[j], 
                                               relativeScale, stepSizeCutoff,
                                               verbose)
                hess[j][i] = hess[i][j]

        return hess


    def ComputeHessianElement(self, costFunc, chiSq, 
                              params, i, j, epsi, epsj, 
                              relativeScale, stepSizeCutoff, verbose):
        return 0.5 * self.hessian_elem(costFunc, chiSq, 
                                       params, i, j, epsi, epsj, 
                                       relativeScale, stepSizeCutoff, 
                                       verbose)

    def CalcHessianInLogParameters(self, params, eps, relativeScale = False, 
                                   stepSizeCutoff = 1e-6, verbose = False):
        return self.hessian_log_params(params, eps, relativeScale,
                                       stepSizeCutoff, verbose)

    def CalcHessian(self, params, epsf, relativeScale = True, 
                    stepSizeCutoff = None, jacobian = None, verbose = False):
    	"""
	Finite difference the residual dictionary to get a dictionary
	for the Hessian. It will be indexed the same as the residuals.
	Note: epsf is either a scalar or an array.
        If relativeScale is False then epsf is the stepsize used (it should
        already be multiplied by typicalValues before Jacobian is called)
        If relativeScale is True then epsf is multiplied by params.
        The two previous statements hold for both scalar and vector valued
        epsf.
        """
        return self.hessian(params, epsf, relativeScale,
                            stepSizeCutoff, jacobian, verbose)

    def CalcHessianUsingResiduals(self,params,epsf,relativeScale = True, moreAcc = False) :
    	currParams = copy.copy(params)
	nOp = len(currParams)
	J,jtj = self.GetJandJtJ(currParams)
	self.CalculateForAllDataPoints(currParams)
	self.ComputeInternalVariables()
	localCalcVals = self.calcVals
	localIntVars = self.internalVars

	paramlist = scipy.array(currParams)

        # epsf may be an array or scalar
        if scipy.isscalar(epsf):
            min_epsf = epsf
        else:
            min_epsf = min(epsf)

	if relativeScale is True :
            ## eps should be epsf*abs(typicalvalue) JJW 7.22.05
            eps = epsf * abs(paramlist)
            for i in range(0,len(eps)) :
                if eps[i] < min_epsf:
                    eps[i] = min_epsf
	else:
            eps = epsf * scipy.ones(len(paramlist),scipy.Float)

	secondDeriv = scipy.zeros((nOp,nOp),scipy.Float)
	
	if moreAcc == False :
		for index in range(0,len(params)) :
            		paramsPlus = currParams.__copy__()
	    		paramsPlus[index] = paramsPlus[index] + eps[index]
	    		JPlus = self.GetJacobian(paramsPlus)
	    		for res in J.keys() :
				resvalue = self.residuals.getByKey(res).GetValue\
					(localCalcVals, localIntVars, currParams)
				secondDeriv[index,:] += (scipy.array(JPlus[res])-scipy.array(J[res]))/ \
					eps[index]*resvalue

	elif moreAcc == True :
		for index in range(0,len(params)) :
                        paramsPlus = currParams.__copy__()
                        paramsPlus[index] = paramsPlus[index] + eps[index]
                        JPlus = self.GetJacobian(paramsPlus)
                       	paramsMinus = currParams.__copy__()
			paramsMinus[index] = paramsMinus[index] - eps[index] 
			JMinus =  self.GetJacobian(paramsMinus)
	
			for res in J.keys() :
                                resvalue = self.residuals.getByKey(res).GetValue\
                                        (localCalcVals, localIntVars, currParams)
                                secondDeriv[index,:] += (scipy.array(JPlus[res])-scipy.array(JMinus[res]))/ \
                                        (2.0*eps[index])*resvalue
	
	self.params.update(currParams)
	# force symmetry, because it will never be exactly symmetric.
	# Should check that (jtj+secondDeriv) is approximately symmetric if concerned
	# about accuracy of the calculation --- the difference between the i,j and the
	# j,i entries should give you an indication of the accuracy
	hess = 0.5*(scipy.transpose(jtj + secondDeriv) + (jtj + secondDeriv))
	return hess

    def CalcHessianUsingResidualsInLogParams(self,params,epsf,relativeScale = True, moreAcc = False) :
    	currParams = copy.copy(params)
        currRawParams = copy.copy(params)
        currRawParams.update(scipy.exp(params))
	nOp = len(currParams)
	J,jtj = self.GetJandJtJInLogParameters(currParams)
	self.ComputeInternalVariables()
	localCalcVals = self.calcVals
	localIntVars = self.internalVars

	paramlist = scipy.array(currParams)
        eps = epsf * scipy.ones(len(paramlist),scipy.Float)

	secondDeriv = scipy.zeros((nOp,nOp),scipy.Float)
	if moreAcc == False :	
		for index in range(0,len(params)) :
            		paramsPlus = currParams.__copy__()
	    		paramsPlus[index] = paramsPlus[index] + eps[index]
	    		JPlus,tmp = self.GetJandJtJInLogParameters(paramsPlus)

	    		for res in J.keys() :
				resvalue = self.residuals.getByKey(res).GetValue\
					(localCalcVals, localIntVars, currRawParams)
				secondDeriv[index,:] += (scipy.array(JPlus[res])-scipy.array(J[res]))/ \
					eps[index]*resvalue

	elif moreAcc == True :
		for index in range(0,len(params)) :
                        paramsPlus = currParams.__copy__()
                        paramsPlus[index] = paramsPlus[index] + eps[index]
                        JPlus,tmp = self.GetJandJtJInLogParameters(paramsPlus)
			paramsMinus = currParams.__copy__()
                        paramsMinus[index] = paramsMinus[index] - eps[index]
                       	JMinus,tmp = self.GetJandJtJInLogParameters(paramsMinus) 
			
			for res in J.keys() :
                                resvalue = self.residuals.getByKey(res).GetValue\
                                        (localCalcVals, localIntVars, currRawParams)
                                print JPlus[res], '\n\n', JMinus[res], '\n\n', eps[index], '\n\n', res, resvalue 
                                secondDeriv[index,:] += (scipy.array(JPlus[res])-scipy.array(JMinus[res]))/ \
                                        (2.0*eps[index])*resvalue
	self.params.update(scipy.exp(currParams))
	# force symmetry, because it will never be exactly symmetric.
	# Should check that (jtj+secondDeriv) is approximately symmetric if concerned
	# about accuracy of the calculation --- the difference between the i,j and the
	# j,i entries should give you an indication of the accuracy
	hess = 0.5*(scipy.transpose(jtj + secondDeriv) + (jtj + secondDeriv))
	return hess
        
    def CalcResidualResponseArray(self, j, h) :
        """
        Calculate the Residual Response array. This array represents the change in
        a residual obtained by a finite change in a data value.

        Inputs:
        (self, j, h)
        j -- jacobian matrix to use
        h -- hessian matrix to use

        Outputs:
        response -- The response array
        """
        j,h = scipy.asarray(j), scipy.asarray(h)
	[m,n] = j.shape
	response = scipy.zeros((m,m),scipy.Float)
	ident = scipy.eye(m,typecode=scipy.Float)
	hinv = scipy.linalg.pinv2(h,1e-40)
	tmp = scipy.matrixmultiply(hinv,scipy.transpose(j))
	tmp2 = scipy.matrixmultiply(j,tmp)
	response = ident - tmp2

	return response

    def CalcParameterResponseToResidualArray(self,j,h):
        """
        Calculate the parameter response to residual array. This array represents
        the change in parameter resulting from a change in data (residual).

        Inputs:
        (self, j, h)
        j -- jacobian matrix to use
        h -- hessian matrix to use

        Outputs:
        response -- The response array
        """
        j,h = scipy.asarray(j), scipy.asarray(h)
	[m,n] = j.shape
	response = scipy.zeros((n,m),scipy.Float)
	hinv = scipy.linalg.pinv2(h,1e-40)
	response = -scipy.matrixmultiply(hinv,scipy.transpose(j))

	return response
        
    ############################################################################
    # Getting/Setting variables below

    def SetExperimentCollection(self, exptColl):
        self.exptColl = exptColl

        for exptKey, expt in exptColl.items():
            exptData = expt.GetData()
            for calcKey, calcData in exptData.items():
                for depVarKey, depVarData in calcData.items():
                    sortedData = depVarData.items()
                    sortedData.sort()
                    for indVar, (value, uncert) in sortedData:
                        resName = (exptKey, calcKey, depVarKey, indVar)
                        res = Residuals.ScaledErrorInFit(resName, depVarKey,
                                                         calcKey, indVar, value,
                                                         uncert, exptKey)
                        self.residuals.setByKey(resName, res)

            # Add in the PeriodChecks
            for period in expt.GetPeriodChecks():
                calcKey, depVarKey, indVarValue = period['calcKey'], period['depVarKey'], period['startTime']
                resName = (exptKey, calcKey, depVarKey, indVarValue, 'PeriodCheck')
                res = Residuals.PeriodCheckResidual(resName, calcKey, depVarKey, indVarValue,
                                                    period['period'], period['sigma'])
                self.residuals.setByKey(resName, res)

            # Add in the AmplitudeChecks
            for amplitude in expt.GetAmplitudeChecks():
                calcKey, depVarKey = amplitude['calcKey'], amplitude['depVarKey']
                indVarValue0, indVarValue1 = amplitude['startTime'], amplitude['testTime']
                resName = (exptKey, calcKey, depVarKey, indVarValue0, indVarValue1, 'AmplitudeCheck')
                res = Residuals.AmplitudeCheckResidual(resName, calcKey, depVarKey, indVarValue0, indVarValue1,
                                                       amplitude['period'], amplitude['sigma'], exptKey)
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

def end_pypar():
    """
    Used by atexit to call pypar.finalize only at the end,
    and only once.
    """
    if HAVE_PYPAR:
        import pypar
        pypar.finalize()

import atexit
atexit.register(end_pypar)
