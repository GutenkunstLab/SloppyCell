"""
Model class that unites theory with data.
"""

import logging
logger = logging.getLogger('Model_mod')

import copy
import sets

import scipy

import SloppyCell
import SloppyCell.Residuals as Residuals
import SloppyCell.Collections as Collections
import SloppyCell.Utility as Utility
import KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

_double_epsilon_ = scipy.finfo(scipy.float_).eps
    
class Model:
    """
    A Model object connects a set of experimental data with the objects used to
    model that data.

    Most importantly, a Model can calculate a cost for a given set of 
    parameters, characterizing how well those parameters fit the data contained
    within the model.
    """

    imag_cutoff = 1e-8

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
        elif isinstance(calcs, dict):
            calcs = Collections.CalculationCollection(calcs.values())
        self.SetCalculationCollection(calcs)

        self.observers = KeyedList()
        self.parameter_bounds = {}

    def compile(self):
        """
        Compile all the calculations contained within the Model.
        """
        for calc in self.get_calcs().values():
            calc.compile()

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
            self.calcColl.get(calcName).set_var_ic(varName, initialValue)

    def _evaluate(self, params, T=1):
        """
        Evaluate the cost for the model, returning the intermediate residuals,
        and chi-squared.

        (Summing up the residuals is a negligible amount of work. This 
         arrangment makes notification of observers much simpler.)
        """
        self.params.update(params)
        self.check_parameter_bounds(params)
        self.CalculateForAllDataPoints(params)
        self.ComputeInternalVariables(T)

        resvals = [res.GetValue(self.calcVals, self.internalVars, self.params)
                   for res in self.residuals.values()]

        # Occasionally it's useful to use residuals with a sqrt(-1) in them,
        #  to get negative squares. Then, however, we might get small imaginary
        #  parts in our results, which this shaves off.
        chisq = scipy.real_if_close(scipy.sum(scipy.asarray(resvals)**2), 
                                    tol=self.imag_cutoff)
        if scipy.isnan(chisq):
            logger.warn('Chi^2 is NaN, converting to Infinity.')
            chisq = scipy.inf
        cost = 0.5 * chisq

        entropy = 0
        for expt, sf_ents in self.internalVars['scaleFactor_entropies'].items():
            for group, ent in sf_ents.items():
                entropy += ent

        self._notify(event = 'evaluation', 
                     resvals = resvals,
                     chisq = chisq,
                     cost = cost, 
                     free_energy = cost-T*entropy,
                     entropy = entropy,
                     params = self.params)

        return resvals, chisq, cost, entropy

    def res(self, params):
        """
        Return the residual values of the model fit given a set of parameters
        """

        return self._evaluate(params)[0]

    def res_log_params(self, log_params):
        """
        Return the residual values given the logarithm of the parameters
        """
        return self.res(scipy.exp(log_params))

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

    def cost_log_params(self, log_params):
        """
        Return the cost given the logarithm of the input parameters
        """
        return self.cost(scipy.exp(log_params))

    def free_energy(self, params, T):
        temp, temp, c, entropy = self._evaluate(params, T=T)
        return c - T * entropy

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
	
    resDict = res_dict
    # ...


    def AddResidual(self, res):
        self.residuals.setByKey(res.key, res)

    def Force(self, params, epsf, relativeScale=False, stepSizeCutoff=None):
        """
        Force(parameters, epsilon factor) -> list

        Returns a list containing the numerical gradient of the cost with
        respect to each parameter (in the parameter order of the
        CalculationCollection). Each element of the gradient is:
            cost(param + eps) - cost(param - eps)/(2 * eps).
        If relativeScale is False then epsf is the stepsize used (it should
        already be multiplied by typicalValues before Jacobian is called)
        If relativeScale is True then epsf is multiplied by params.
        The two previous statements hold for both scalar and vector valued
        epsf.
        """

        force = []
        params = scipy.array(params)
        
        if stepSizeCutoff==None:
            stepSizeCutoff = scipy.sqrt(_double_epsilon_)
            
	if relativeScale is True:
            eps = epsf * abs(params)
	else:
            eps = epsf * scipy.ones(len(params),scipy.float_)

        for i in range(0,len(eps)):
            if eps[i] < stepSizeCutoff:
                eps[i] = stepSizeCutoff

        for index, param in enumerate(params):
            paramsPlus = params.copy()
            paramsPlus[index] = param + eps[index]
            costPlus = self.cost(paramsPlus)

            paramsMinus = params.copy()
            paramsMinus[index] = param - eps[index]
            costMinus = self.cost(paramsMinus)

            force.append((costPlus-costMinus)/(2.0*eps[index]))

        return force

    def gradient_sens(self, params):
        """
        Return the gradient of the cost, d_cost/d_param as a KeyedList.

        This method uses sensitivity integration, so it only applies to
        ReactionNetworks.
        """
        self.params.update(params)

        # The cost is 0.5 * sum(res**2), 
        # so the gradient is sum(res * dres_dp)

        jac_dict = self.jacobian_sens(params)
        res_dict = self.res_dict(params)

        force = scipy.zeros(len(params), scipy.float_)
        for res_key, res_val in res_dict.items():
            res_derivs = jac_dict.get(res_key)
            force += res_val * scipy.asarray(res_derivs)

        gradient = self.params.copy()
        gradient.update(force)

        return gradient

    def gradient_log_params_sens(self, log_params):
        """
        Return the gradient of the cost wrt log parameters, d_cost/d_log_param
        as a KeyedList.

        This method uses sensitivity integration, so it only applies to
        ReactionNetworks.
        """
        # We just need to multiply dcost_dp by p.
        params = scipy.exp(log_params)
        gradient = self.gradient_sens(params)
        gradient_log = gradient.copy()
        gradient_log.update(scipy.asarray(gradient) * scipy.asarray(params))

        return gradient_log

    def CalculateForAllDataPoints(self, params):
        """
        CalculateForAllDataPoints(parameters) -> dictionary

        Gets a dictionary of measured independent variables indexed by 
        calculation from the ExperimentCollection and passes it to the
        CalculationCollection. The returned dictionary is of the form:
         dictionary[experiment][calculation][dependent variable]
                   [independent variabled] -> calculated value.
        """
        self.params.update(params)

        varsByCalc = self.GetExperimentCollection().GetVarsByCalc()
        self.calcVals = self.GetCalculationCollection().Calculate(varsByCalc, 
                                                                  params)
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
        self.calcVals, self.calcSensitivityVals =\
                self.GetCalculationCollection().CalculateSensitivity(varsByCalc,
                                                                     params)
        return self.calcSensitivityVals

    def ComputeInternalVariables(self, T=1):
        sf, sf_ents = self.compute_scale_factors(T)
        self.internalVars['scaleFactors'] = sf
        self.internalVars['scaleFactor_entropies'] = sf_ents

    def compute_scale_factors(self, T):
        """
        Compute the scale factors for the current parameters and return a dict.

        The dictionary is of the form dict[exptId][varId] = scale_factor
        """
        scale_factors = {}
        scale_factor_entropies = {}
        for exptId, expt in self.GetExperimentCollection().items():
            scale_factors[exptId], scale_factor_entropies[exptId] =\
                    self._compute_sf_and_sfent_for_expt(expt, T)

        return scale_factors, scale_factor_entropies

    def _compute_sf_and_sfent_for_expt(self, expt, T):
        # Compute the scale factors for a given experiment
        scale_factors = {}
        scale_factor_entropies = {}

        exptData = expt.GetData()
        expt_integral_data = expt.GetIntegralDataSets()
        fixed_sf = expt.get_fixed_sf()

        sf_groups = expt.get_sf_groups()

        for group in sf_groups:
            # Do any of the variables in this group have fixed scale factors?
            fixed = sets.Set(group).intersection(sets.Set(fixed_sf.keys()))
            fixedAt = sets.Set([fixed_sf[var] for var in fixed])

            # We'll need to index the scale factor entropies on the *group*
            # that shares a scale factor, since we only have one entropy per
            # shared scale factor. So we need to index on the group of
            # variables. We sort the group and make it hashable to avoid any
            # double-counting.
            hash_group = expt._hashable_group(group)
            if len(fixedAt) == 1:
                value = fixedAt.pop()
                for var in group:
                    scale_factors[var] = value
                    scale_factor_entropies[hash_group] = 0
                continue
            elif len(fixedAt) > 1:
                    raise ValueError('Shared scale factors fixed at '
                                     'inconsistent values in experiment '
                                     '%s!' % expt.GetName())

            # Finally, compute the scale factor for this group
            theoryDotData, theoryDotTheory = 0, 0
            # For discrete data
            for calc in exptData:
                # Pull out the vars we have measured for this calculation
                for var in sets.Set(group).intersection(sets.Set(exptData[calc].keys())):
                    for indVar, (data, error) in exptData[calc][var].items():
                        theory = self.calcVals[calc][var][indVar]
                        theoryDotData += (theory * data) / error**2
                        theoryDotTheory += theory**2 / error**2

            # Now for integral data
            for dataset in expt_integral_data:
                calc = dataset['calcKey']
                theory_traj = self.calcVals[calc]['full trajectory']
                data_traj = dataset['trajectory']
                uncert_traj = dataset['uncert_traj']
                interval = dataset['interval']
                T = interval[1] - interval[0]
                for var in group.intersection(sets.Set(dataset['vars'])):
                    TheorDotT = self._integral_theorytheory(var, theory_traj,
                                                            uncert_traj, 
                                                            interval)
                    theoryDotTheory += TheorDotT/T
                    TheorDotD = self._integral_theorydata(var, theory_traj,
                                                          data_traj, 
                                                          uncert_traj, 
                                                          interval)
                    theoryDotData += TheorDotD/T

            # Now for the extrema data
            for ds in expt.scaled_extrema_data:
                calc = ds['calcKey']
                if ds['type']  == 'max':
                    var = ds['var'] + '_maximum'
                elif ds['type'] == 'min':
                    var = ds['var'] + '_minimum'
                data, error = ds['val'], ds['sigma']
                theory = self.calcVals[calc][var]\
                        [ds['minTime'],ds['maxTime']][1]
                theoryDotData += (theory * data) / error**2
                theoryDotTheory += theory**2 / error**2

            for var in group:
                if theoryDotTheory != 0:
                    scale_factors[var] = theoryDotData/theoryDotTheory
                else:
                    scale_factors[var] = 1
                entropy = expt.compute_sf_entropy(hash_group, theoryDotTheory,
                                                  theoryDotData, T)
                scale_factor_entropies[hash_group] = entropy

        return scale_factors, scale_factor_entropies

    def _integral_theorytheory(self, var, theory_traj, uncert_traj, interval):
        def integrand(t):
            theory = theory_traj.evaluate_interpolated_traj(var, t)
            uncert = uncert_traj.evaluate_interpolated_traj(var, t)
            return theory**2/uncert**2
        val, error = scipy.integrate.quad(integrand, interval[0], interval[1],
                                          limit=int(1e5))
        return val

    def _integral_theorydata(self, var, theory_traj, data_traj, uncert_traj, 
                             interval):
        def integrand(t):
            theory = theory_traj.evaluate_interpolated_traj(var, t)
            data = data_traj.evaluate_interpolated_traj(var, t)
            uncert = uncert_traj.evaluate_interpolated_traj(var, t)
            return theory*data/uncert**2
        val, error = scipy.integrate.quad(integrand, interval[0], interval[1],
                                          limit=int(1e5))
        return val

    def ComputeInternalVariableDerivs(self):
        """
        compute_scale_factorsDerivs() -> dictionary

        Returns the scale factor derivatives w.r.t. parameters
	appropriate for each chemical in each
        experiment, given the current parameters. The returned dictionary
        is of the form: internalVarsDerivs['scaleFactors'] \
                = dict[experiment][chemical][parametername] -> derivative.
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
            # Now for the extrema data
            for ds in expt.scaled_extrema_data:
                exptDepVars.add(ds['var'])

            for depVar in exptDepVars:
		self.internalVarsDerivs['scaleFactors'][exptName][depVar] = {}
		if depVar in expt.GetFixedScaleFactors():
                    for pname in p.keys():
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

                for ds in expt.scaled_extrema_data:
                    if ds['type'] == 'max':
                        var = ds['var'] + '_maximum'
                    elif ds['type'] == 'min':
                        var = ds['var'] + '_minimum'
                    data, error = ds['val'], ds['sigma']
                    theory = self.calcVals[ds['calcKey']][var]\
                            [ds['minTime'],ds['maxTime']][1]
                    theoryDotData += (theory * data) / error**2
                    theoryDotTheory += theory**2 / error**2

                # now get derivative of the scalefactor
                for pname in p.keys():
                    theorysensDotData, theorysensDotTheory = 0, 0

                    for calc in exptData:
                        clc = self.calcColl.get(calc)
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
                    for ds in expt.scaled_extrema_data:
                        if ds['type'] == 'max':
                            var = ds['var'] + '_maximum'
                        elif ds['type'] == 'min':
                            var = ds['var'] + '_minimum'
                        theory = self.calcVals[ds['calcKey']][var]\
                                [ds['minTime'],ds['maxTime']][1]
                        data, error = ds['val'], ds['sigma']
                        theorysens = self.calcSensitivityVals[ds['calcKey']][var][ds['minTime'],ds['maxTime']].get(pname, 0.0)
                        theorysensDotData += (theorysens * data) / error**2
                        theorysensDotTheory += (theorysens * theory) / error**2

                    deriv_dict = self.internalVarsDerivs['scaleFactors'][exptName][depVar]
                    try:
                        deriv_dict[pname] = theorysensDotData/theoryDotTheory\
                                - 2*theoryDotData*theorysensDotTheory/(theoryDotTheory)**2
                    except ZeroDivisionError:
                        deriv_dict[pname] = 0

	return self.internalVarsDerivs['scaleFactors']

    def jacobian_log_params_sens(self, log_params):
    	"""
	Return a KeyedList of the derivatives of the model residuals w.r.t.
        the lograithms of the parameters parameters.

        The method uses the sensitivity integration. As such, it will only
        work with ReactionNetworks.

        The KeyedList is of the form:
            kl.get(resId) = [dres/dlogp1, dres/dlogp2...]
	"""
        params = scipy.exp(log_params)
        j = self.jacobian_sens(params)
        j_log = j.copy()
        j_log.update(scipy.asarray(j) * scipy.asarray(params))

        return j_log

    def jacobian_sens(self, params):
    	"""
	Return a KeyedList of the derivatives of the model residuals w.r.t.
        parameters.

        The method uses the sensitivity integration. As such, it will only
        work with ReactionNetworks.

        The KeyedList is of the form:
            kl[resId] = [dres/dp1, dres/dp2...]
	"""
        self.params.update(params)

        # Calculate sensitivities
	self.CalculateSensitivitiesForAllDataPoints(params)
	self.ComputeInternalVariables()
	self.ComputeInternalVariableDerivs()

        # Calculate residual derivatives
        deriv = [(resId, res.Dp(self.calcVals, self.calcSensitivityVals,
                                self.internalVars, self.internalVarsDerivs,
                                self.params))
                 for (resId, res) in self.residuals.items()]

	return KeyedList(deriv)

    def jacobian_fd(self, params, eps, 
                    relativeScale=False, stepSizeCutoff=None):
    	"""
	Return a KeyedList of the derivatives of the model residuals w.r.t.
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
            stepSizeCutoff = scipy.sqrt(_double_epsilon_)
            
	if relativeScale:
            eps_l = scipy.maximum(eps * abs(params), stepSizeCutoff)
	else:
            eps_l = scipy.maximum(eps * scipy.ones(len(params),scipy.float_),
                                  stepSizeCutoff)

	J = KeyedList() # will hold the result
	for resId in res.keys():
            J.set(resId, [])
        # Two-sided finite difference
	for ii in range(len(params)):
            params[ii] = orig_vals[ii] + eps_l[ii]
	    resPlus = self.resDict(params)

            params[ii] = orig_vals[ii] - eps_l[ii]
            resMinus = self.resDict(params)

            params[ii] = orig_vals[ii]

	    for resId in res.keys():
                res_deriv = (resPlus[resId]-resMinus[resId])/(2.*eps_l[ii])
                J.get(resId).append(res_deriv)
	
	# NOTE: after call to ComputeResidualsWithScaleFactors the Model's
	# parameters get updated, must reset this:
        self.params.update(params)
	return J

    def GetJacobian(self,params):
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

    def GetJandJtJ(self,params):
	j = self.GetJacobian(params)
	mn = scipy.zeros((len(params),len(params)),scipy.float_)

	for paramind in range(0,len(params)):
	  for paramind1 in range(0,len(params)):
	    sum = 0.0
	    for kys in j.keys():
	      sum = sum + j.get(kys)[paramind]*j.get(kys)[paramind1]

	    mn[paramind][paramind1] = sum
    	return j,mn
   
    def GetJandJtJInLogParameters(self,params):
    	# Formula below is exact if you have perfect data. If you don't
	# have perfect data (residuals != 0) you get an extra term when you
	# compute d^2(cost)/(dlogp[i]dlogp[j]) which is
	# sum_resname (residual[resname] * jac[resname][j] * delta_jk * p[k])
	# but can be ignored when residuals are zeros, and maybe should be
	# ignored altogether because it can make the Hessian approximation 
        # non-positive definite
	pnolog = scipy.exp(params)	
	jac, jtj = self.GetJandJtJ(pnolog)
	for i in range(len(params)):
            for j in range(len(params)):
                jtj[i][j] = jtj[i][j]*pnolog[i]*pnolog[j]

	res = self.resDict(pnolog)	
	for resname in self.residuals.keys():
            for j in range(len(params)):
               	# extra term --- not including it 
		# jtj[j][j] += res[resname]*jac[resname][j]*pnolog[j]
                jac.get(resname)[j] = jac.get(resname)[j]*pnolog[j]

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
            stepSizeCutoff = scipy.sqrt(_double_epsilon_)
            
        params = scipy.asarray(params)
	if relativeScale:
            eps = epsf * abs(params)
	else:
            eps = epsf * scipy.ones(len(params),scipy.float_)

        # Make sure we don't take steps smaller than stepSizeCutoff
        eps = scipy.maximum(eps, stepSizeCutoff)

        if jacobian is not None:
            # Turn off the relative scaling since that would overwrite all this
            relativeScale = False

            jacobian = scipy.asarray(jacobian)
            if len(jacobian.shape) == 0:
                resDict = self.resDict(params)
                new_jacobian = scipy.zeros(len(params),scipy.float_)
                for key, value in resDict.items():
                    new_jacobian += 2.0*value*scipy.array(jacobian[0][key])
                jacobian = new_jacobian
            elif len(jacobian.shape) == 2: # Need to sum up the total jacobian
                residuals = scipy.asarray(self.res(params))
                # Changed by rng7. I'm not sure what is meant by "sum up the
                # total jacobian". The following line failed due to shape
                # mismatch. From the context below, it seems that the dot
                # product is appropriate.
                #jacobian = 2.0*residuals*jacobian
                jacobian = 2.0 * scipy.dot(residuals, jacobian)

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

	hess = scipy.zeros((nOv, nOv), scipy.float_)

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
            eps = scipy.ones(len(params), scipy.float_) * eps

	## compute cost at f(x)
	f0 = self.cost_log_params(scipy.log(params))

	hess = scipy.zeros((nOv, nOv), scipy.float_)

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

    def CalcResidualResponseArray(self, j, h):
        """
        Calculate the Residual Response array. This array represents the change
        in a residual obtained by a finite change in a data value.

        Inputs:
        (self, j, h)
        j -- jacobian matrix to use
        h -- hessian matrix to use

        Outputs:
        response -- The response array
        """
        j,h = scipy.asarray(j), scipy.asarray(h)
	[m,n] = j.shape
	response = scipy.zeros((m,m),scipy.float_)
	ident = scipy.eye(m,typecode=scipy.float_)
	hinv = scipy.linalg.pinv2(h,1e-40)
	tmp = scipy.dot(hinv,scipy.transpose(j))
	tmp2 = scipy.dot(j,tmp)
	response = ident - tmp2

	return response

    def CalcParameterResponseToResidualArray(self,j,h):
        """
        Calculate the parameter response to residual array. This array
        represents the change in parameter resulting from a change in data
        (residual).

        Inputs:
        (self, j, h)
        j -- jacobian matrix to use
        h -- hessian matrix to use

        Outputs:
        response -- The response array
        """
        j,h = scipy.asarray(j), scipy.asarray(h)
	[m,n] = j.shape
	response = scipy.zeros((n,m),scipy.float_)
	hinv = scipy.linalg.pinv2(h,1e-40)
	response = -scipy.dot(hinv,scipy.transpose(j))

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
                calcKey, depVarKey, indVarValue = period['calcKey'], \
                        period['depVarKey'], period['startTime']
                resName = (exptKey, calcKey, depVarKey, indVarValue, 
                           'PeriodCheck')
                res = Residuals.PeriodCheckResidual(resName, calcKey, depVarKey,
                                                    indVarValue, 
                                                    period['period'], 
                                                    period['sigma'])
                self.residuals.setByKey(resName, res)

            # Add in the AmplitudeChecks
            for amplitude in expt.GetAmplitudeChecks():
                calcKey, depVarKey = amplitude['calcKey'], \
                        amplitude['depVarKey']
                indVarValue0, indVarValue1 = amplitude['startTime'],\
                        amplitude['testTime']
                resName = (exptKey, calcKey, depVarKey, indVarValue0, 
                           indVarValue1, 'AmplitudeCheck')
                res = Residuals.AmplitudeCheckResidual(resName, calcKey, 
                                                       depVarKey, indVarValue0,
                                                       indVarValue1,
                                                       amplitude['period'], 
                                                       amplitude['sigma'], 
                                                       exptKey)
                self.residuals.setByKey(resName, res)

            # Add in the integral data
            for ds in expt.GetIntegralDataSets():
                for var in ds['vars']:
                    resName = (exptKey, ds['calcKey'], var, 'integral data')
                    res = Residuals.IntegralDataResidual(resName, var,
                                                          exptKey,
                                                          ds['calcKey'],
                                                          ds['trajectory'],
                                                          ds['uncert_traj'],
                                                          ds['interval'])
                    self.residuals.setByKey(resName, res)

            for ds in expt.scaled_extrema_data:
                ds['exptKey'] = expt.name
                ds['key'] = '%s_%simum_%s_%s' % (ds['var'], ds['type'],
                                                 str(ds['minTime']),
                                                 str(ds['maxTime']))
                res = Residuals.ScaledExtremum(**ds)
                self.AddResidual(res)

               
    def get_expts(self):
        return self.exptColl

    def set_var_optimizable(self, var, is_optimizable):
        for calc in self.get_calcs().values():
            try:
                calc.set_var_optimizable(var, is_optimizable)
            except KeyError:
                pass
        self.params = self.calcColl.GetParameters()

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

    def add_parameter_bounds(self, param_id, pmin, pmax):
        """
        Add bounds on a specific parameter.
        
        Cost evaluations will raise an exception if these bounds are violated.
        """
        self.parameter_bounds[param_id] = pmin, pmax

    def check_parameter_bounds(self, params):
        self.params.update(params)
        for id, (pmin, pmax) in self.parameter_bounds.items():
            if not pmin <= self.params.get(id) <= pmax:
                err = 'Parameter %s has value %f, which is outside of given bounds %f to %f.' % (id, self.params.get(id), pmin, pmax)
                raise Utility.SloppyCellException(err)
