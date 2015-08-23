import logging
logger = logging.getLogger('Collections')

import sets, copy

import scipy

import SloppyCell
import SloppyCell.Utility as Utility
import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList

if SloppyCell.HAVE_PYPAR:
    import pypar

class ExperimentCollection(dict):
    """
    ExperimentCollection(experiment name list)

    An ExperimentCollection unites a collection of experiments. For now, it's
    most important function is to collect group the independent variables by
    Calculation to avoid wasting computer effort redoing calculations.

    Individual experiments can be accessed via dictionary-type indexing.
    """
    def __init__(self, exptList = []):
        for expt in exptList:
            self.AddExperiment(expt)

    def AddExperiment(self, expt):
        """
        LoadExperiments(experiment name list)

        Adds the experiments in the list to the collection
        """
        if expt.GetName() in self:
            raise ValueError("Experiment already has name %s" 
                             % str(expt.GetName()))

        self[expt.GetName()] = expt

    def GetVarsByCalc(self):
        """
        GetIndVarsByCalc() -> dictionary

        Returns a dictionary of all the dependent and independent variables for 
        all the calculations required to compare with the data in all the
        experiments. The dictionary is of the form: 
         dictionary(calculation name) -> ordered list of unique independent
                                         variables
        """
        varsByCalc = {}
        for expt in self.values():
            data = expt.GetData()

            for calc in data:
                varsByCalc.setdefault(calc, {})
                for depVar in data[calc]:
                    # Using a set is a convenient way to make sure
                    # independent variables aren't repeated
                    varsByCalc[calc].setdefault(depVar, sets.Set())
                    varsByCalc[calc][depVar].\
                            union_update(sets.Set(data[calc][depVar].keys()))

            for period in expt.GetPeriodChecks():
                calc, depVar = period['calcKey'], period['depVarKey']
                start, period = period['startTime'], period['period']
                if calc not in varsByCalc.keys():
                    varsByCalc[calc].setdefault(calc, {})
                if depVar not in varsByCalc[calc]:
                    varsByCalc[calc].setdefault(depVar, sets.Set())
                varsByCalc[calc][depVar].union_update([start, start+2.0*period])

            for amplitude in expt.GetAmplitudeChecks():
                calc, depVar = amplitude['calcKey'], amplitude['depVarKey']
                start, test, period = (amplitude['startTime'], 
                                       amplitude['testTime'], 
                                       amplitude['period'])
                if calc not in varsByCalc.keys():
                    varsByCalc[calc].setdefault(calc, {})
                if depVar not in varsByCalc[calc]:
                    varsByCalc[calc].setdefault(depVar, sets.Set())
                varsByCalc[calc][depVar].union_update([start, start+period,
                                                       test, test+period])

            for data_set in expt.GetIntegralDataSets():
                calc = data_set['calcKey']
                varsByCalc.setdefault(data_set['calcKey'], {})
                varsByCalc[calc][('full trajectory')] = data_set['interval']

            for ds in expt.scaled_extrema_data:
                calc = ds['calcKey']
                varsByCalc.setdefault(ds['calcKey'], {})
                if ds['type'] == 'max':
                    called = ds['var'] + '_maximum'
                elif ds['type'] == 'min': 
                    called = ds['var'] + '_minimum'
                varsByCalc[calc][called] = (ds['minTime'], ds['maxTime'])
                
        # But I convert the sets back to sorted lists before returning
        for calc in varsByCalc:
            for depVar in varsByCalc[calc]:
                varsByCalc[calc][depVar] = list(varsByCalc[calc][depVar])
                varsByCalc[calc][depVar].sort()

        return varsByCalc

    def GetData(self):
        """
        GetData() -> dictionary

        Returns a dictionary containing all the data for the experiments. The
        dictionary is of the form:
         dictionary[expt name][calc name][dependent vars][independent vars]
                 = value.

        Note that value may be an arbitrary object.
        """

        data = {}
        for exptName, expt in self.items():
            data[exptName] = expt.GetData()
        
        return data

class Experiment:
    def __init__(self, name = '', data = {}, fixedScaleFactors = {},
                 longName = '', shared_sf = []):
        self.SetName(name)
        self.SetData(data)
        self.SetFixedScaleFactors(fixedScaleFactors)
        self.set_shared_sf(shared_sf)
        self.periodChecks=[]
        self.amplitudeChecks=[]
        self.integral_data=[]
        self.sf_priors = {}
        self.scaled_extrema_data = []

    def _hashable_group(self, group):
        """
        Return a sorted tuple of the elements of group.
        """
        if isinstance(group, str):
            group = [group]
        hash_group = sets.Set(group)
        hash_group = list(hash_group)
        list(group).sort()
        hash_group = tuple(group)
        return hash_group

    def get_sf_groups(self):
        """
        Return tuples representing all the scale factors in this experiment.

        A tuple will contain multiple entries if several variables share a
        scale factor.
        """
        # Get the variables measured in this experiment
        measuredVars = sets.Set()
        for calcId in self.data:
            measuredVars.union_update(sets.Set(self.data[calcId].keys()))
        for dataset in self.integral_data:
            measuredVars.union_update(sets.Set(dataset['vars']))
        for dataset in self.scaled_extrema_data:
            measuredVars.add(dataset['var'])

        sf_groups = [self._hashable_group(group) for group 
                     in self.get_shared_sf()]
        # Flatten out the list of shared scale factors so we can also get
        #  the unshared ones...
        flattened = []
        for g in sf_groups:
            flattened.extend(g)
        # These are variables that don't share scale factors
        unshared = [self._hashable_group(var) for var in measuredVars
                    if var not in flattened]
        sf_groups.extend(unshared)
        return sf_groups

    def set_sf_prior(self, group, prior_type, prior_params=None):
        """
        Set the type of prior to place on a given group of scalefactors.

        The group contains a collection of variables that are sharing a given
        scale factor which may be just one variable. You can see what the
        current groups are with expt.get_sf_groups().

        Currently implemented prior types are:
            'uniform in sf': This is a uniform prior over scale factors. This
            is simplest and fastest to compute, but it tends to weight
            parameter sets that yield large scale factors heavily. It takes no
            parameters.

            'gaussian in log sf': This is a Gaussian prior over the logarithm
            of the scale factor. This should avoid the problem of weighting
            large factors heavily. It takes two parameters: the mean of the
            normal distribution, and it's standard deviation. For example,
            parameters (log(3.0), log(10)), will place a prior that holds 95%
            of the probability between 3 / 10**2 and 3 * 10**2.
        """
        hash_group = self._hashable_group(group)

        if hash_group not in self.get_sf_groups():
            raise ValueError('Unrecognized group to set scale factor prior on. '
                             'If it is a shared scale factor, you need to '
                             'specify every member of the sharing group.')

        available_types = ['uniform in sf', 'gaussian in log sf']
        if prior_type not in available_types:
            raise ValueError('Unknown prior type %s. Available types are %s.'
                             % (prior_type, available_types))
        self.sf_priors[hash_group] = (prior_type, prior_params)

    def compute_sf_entropy(self, sf_group, theoryDotTheory, theoryDotData, T):
        """
        Compute the entropy for a given scale factor.
        """
        try:
            prior_type, prior_params = self.sf_priors[sf_group]
        except KeyError:
            prior_type, prior_params = 'uniform in sf', None

        if prior_type == 'uniform in sf':
            if theoryDotTheory != 0:
                entropy = scipy.log(scipy.sqrt(2*scipy.pi*T/theoryDotTheory))
            else:
                entropy = 0
        elif prior_type == 'gaussian in log sf':
            # This is the integrand for the prior. Note that it's important
            #  that u = 0 corresponds to B_best. This ensures that the
            #  integration doesn't miss the (possibly sharp) peak there.
            try:
                import SloppyCell.misc_c
                integrand = SloppyCell.misc_c.log_gaussian_prior_integrand
            except ImportError:
                logger.warn('Falling back to python integrand on log gaussian '
                            'prior integration.')
                exp = scipy.exp
                def integrand(u, ak, bk, mulB, siglB, T, B_best, lB_best):
                    B = exp(u) * B_best
                    lB = u + lB_best
                    ret = exp(-ak/(2*T) * (B-B_best)**2
                              - (lB-mulB)**2/(2 * siglB**2))
                    return ret

            mulB, siglB = prior_params
            B_best = theoryDotData/theoryDotTheory
            lB_best = scipy.log(B_best)

            # Often we get overflow errors in the integration, but they
            #  don't actually cause a problem. (exp(-inf) is still 0.) This
            #  prevents scipy from printing annoying error messages in this
            #  case.
            prev_err = scipy.seterr(over='ignore')
            int_args = (theoryDotTheory, theoryDotData, mulB, siglB, T, B_best,
                        lB_best)
            ans, temp = scipy.integrate.quad(integrand, -scipy.inf, scipy.inf, 
                                             args = int_args, limit=1000)
            scipy.seterr(**prev_err)
            entropy = scipy.log(ans)
        else:
            raise ValueError('Unrecognized prior type: %s.' % prior_type)

        return entropy
        
    def SetName(self, name):
        self.name = name

    def GetName(self):
        return self.name

    def set_data(self, data):
        self.data = copy.copy(data)

    def update_data(self, newData):
        self.data.update(newData)

    def get_data(self):
        return self.data

    def set_fixed_sf(self, fixed_sf):
        self.fixedScaleFactors = fixed_sf

    def set_shared_sf(self, shared_sf):
        self.shared_sf = shared_sf

    def get_shared_sf(self):
        return self.shared_sf

    def get_fixed_sf(self):
        return self.fixedScaleFactors


    SetData = set_data
    GetData = get_data
    UpdateData = update_data
    SetFixedScaleFactors = set_fixed_sf
    GetFixedScaleFactors = get_fixed_sf

    def AddPeriodCheck(self, calcKey, chemical, period, sigma, startTime=0.0):
        """
        Constrain the period of the oscillations to a value (period)
        with the error (sigma). The period is found using the maximum
        to maximum distance of the first two maxima found between
        startTime and two periods after the startTime.
        """
        self.periodChecks.append({'calcKey':calcKey, 'depVarKey':chemical,
                                  'period': period, 'sigma': sigma, 
                                  'startTime': startTime})

    def GetPeriodChecks(self):
        return self.periodChecks

    def AddAmplitudeCheck(self, calcKey, chemical, startTime, testTime, period,
                          sigma):
        """
        Turn on applying a constraint that the integrated
        area in two different parts of the plot should be the
        same. startTime and testTime are the starting points to
        begin the integration for the period-long each.
        """
        self.amplitudeChecks.append({'calcKey': calcKey, 'depVarKey': chemical,
                                     'startTime': startTime, 
                                     'testTime': testTime,
                                     'period': period, 'sigma': sigma})

    def GetAmplitudeChecks(self):
        return self.amplitudeChecks

    def GetIntegralDataSets(self):
        return self.integral_data

    def add_integral_data(self, calcKey, traj, uncert_traj, vars, 
                          interval=None):
        """
        Add an integral data set to the experiment.

        calcKey      The id of the calculation this data corresponds to
        traj         The trajectory to compare against
        uncert_traj  A trajectory of data uncertainties
        vars         What variables to fit against
        interval     The time interval to fit over, defaults to the entire traj
        """
        traj.build_interpolated_traj()
        uncert_traj.build_interpolated_traj()
        if interval is None:
            interval = (traj.get_times()[0], traj.get_times()[-1])
        self.integral_data.append({'calcKey': calcKey, 'trajectory': traj,
                                   'uncert_traj': uncert_traj, 'vars': vars,
                                   'interval': interval})

    def add_scaled_max(self, calcKey, var, maxval, sigma, 
                           minTime=None, maxTime=None):
        self.scaled_extrema_data.append({'calcKey': calcKey, 'var':var,
                                         'val':maxval, 'sigma':sigma,
                                         'minTime': minTime, 'maxTime':maxTime,
                                         'type':'max'})
    def add_scaled_min(self, calcKey, var, minval, sigma, 
                           minTime=None, maxTime=None):
        self.scaled_extrema_data.append({'calcKey': calcKey, 'var':var,
                                         'val':minval, 'sigma':sigma,
                                         'minTime': minTime, 'maxTime':maxTime,
                                         'type':'min'})

class CalculationCollection(KeyedList):
    """
    CalculationCollection(calculation name list)

    An CalculationCollection unites a collection of calculations. It is
    responsible for generating a master list of parameters and passing each
    Calculation its appropriate parameters.

    Individual calculations can be accessed via dictionary-type indexing.
    
    XXX: Note that the parameter shuffling has not been extensively tested.
    """

    def __init__(self, calcList = []):
        KeyedList.__init__(self)

        self.params = KeyedList()
        for calc in calcList:
            try:
                if len(calc) == 2:
                    self.AddCalculation(calc[1])
                else:
                    raise ValueError('Incorrect form for calcList')
            except:
                self.AddCalculation(calc)

    def AddCalculation(self, calc):
        """
        LoadCalculations(calculations name list)

        Adds the calculations in the list to the collection and adds their
        parameters to the parameterSet
        """
        if calc.GetName() in self:
            raise ValueError("Calculation already has name %s" 
                             % str(calc.GetName()))

        self.set(calc.GetName(), calc )
        for pName, pValue in calc.GetParameters().items():
            self.params.setdefault(pName, pValue)

    def Calculate(self, varsByCalc, params = None):
        """
        Calculate model predictions for everything in varsByCalc.

        varsByCalc is a dictionary of the form:
            dict[calc name][dep var] = ind var
        
        The return dictionary is of the form:
            dictionary[calc name][dep var][ind var] = result
        """
        if params is not None:
            self.params.update(params)

        results = {}

        calcs_to_do = varsByCalc.keys()
        # Record which calculation each node is doing
        calc_assigned = {}
        while calcs_to_do:
            # The number of calculations to do this round. We want to use
            #  all the processors if possible.
            len_this_block = min(SloppyCell.num_procs, len(calcs_to_do))

            for worker in range(1, len_this_block):
                calc = calcs_to_do.pop()
                calc_assigned[worker] = calc
                logger.debug('Assigning calculation %s to worker %i.'
                             % (calc, worker))
                command = 'Network.calculate(net, vars, params)'
                args = {'net': self.get(calc), 'vars': varsByCalc[calc],
                        'params': self.params}
                pypar.send((command, args), worker)

            # The master does his share here
            calc = calcs_to_do.pop()
            # We use the finally statement because we want to ensure that we
            #  *always* wait for replies from the workers, even if the master
            #  encounters an exception in his evaluation.
            try:
                results[calc] = self.get(calc).calculate(varsByCalc[calc], 
                                                         self.params)
            finally:
                # Collect results from the workers
                for worker in range(1, len_this_block):
                    logger.debug('Receiving result from worker %i.' % worker)
                    results[calc_assigned[worker]] = pypar.receive(worker)
                # If the master encounts an exception, we'll break out of the
                #  function ***here***

            # Check the results we received. If any is a SloppyCellException, 
            #  reraise it.
            for val in results.values():
                if isinstance(val, Utility.SloppyCellException):
                    raise val

        return results

    def CalculateSensitivity(self, varsByCalc, params = None):
        """
        Calculate sensitivities for model predictions of everything in 
        varsByCalc.

        varsByCalc is a dictionary of the form:
            dict[calc name][dep var] = ind var
        
        The return dictionary is of the form:
            dictionary[calc name][dep var][ind var][param] = result
        """
        if params is not None :
            self.params.update(params)
		
        calcSensVals, calcVals = {}, {}
        for (calcName, vars) in varsByCalc.items():
            calc = self.get(calcName)
            vars = varsByCalc[calcName]
            calcPOrder = calc.GetParameters().keys()
            calc.CalculateSensitivity(varsByCalc[calcName], self.params)
            calcSensVals[calcName] = calc.GetSensitivityResult(vars)
            calcVals[calcName] = calc.GetResult(vars)

        return calcVals, calcSensVals

    def GetParameters(self):
        """
        Return a deep copy of the collections parameter KeyedList.
        """
        self.params = KeyedList()
        for calc in self.values():
            for pName, pValue in calc.GetParameters().items():
                self.params.setdefault(pName, pValue)
        return self.params
