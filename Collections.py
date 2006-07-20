import logging
logger = logging.getLogger('Collections')

import sets, copy

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
            raise ValueError, "Experiment already has name %s" % str(expt.GetName())

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
                start, test, period = amplitude['startTime'], amplitude['testTime'], amplitude['period']
                if calc not in varsByCalc.keys():
                    varsByCalc[calc].setdefault(calc, {})
                if depVar not in varsByCalc[calc]:
                    varsByCalc[calc].setdefault(depVar, sets.Set())
                varsByCalc[calc][depVar].union_update([start, start+period,
                                                       test, test+period])
                
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
                                  'period': period, 'sigma': sigma, 'startTime': startTime})

    def GetPeriodChecks(self):
        return self.periodChecks

    def AddAmplitudeCheck(self, calcKey, chemical, startTime, testTime, period, sigma):
        """
        Turn on applying a constraint that the integrated
        area in two different parts of the plot should be the
        same. startTime and testTime are the starting points to
        begin the integration for the period-long each.
        """
        self.amplitudeChecks.append({'calcKey': calcKey, 'depVarKey': chemical,
                                   'startTime': startTime, 'testTime': testTime,
                                   'period': period, 'sigma': sigma})

    def GetAmplitudeChecks(self):
        return self.amplitudeChecks
        

class CalculationCollection(dict):
    """
    CalculationCollection(calculation name list)

    An CalculationCollection unites a collection of calculations. It is
    responsible for generating a master list of parameters and passing each
    Calculation its appropriate parameters.

    Individual calculations can be accessed via dictionary-type indexing.
    
    XXX: Note that the parameter shuffling has not been extensively tested.
    """

    def __init__(self, calcList = []):
        self.params = KeyedList()
        for calc in calcList:
            self.AddCalculation(calc)

    def AddCalculation(self, calc):
        """
        LoadCalculations(calculations name list)

        Adds the calculations in the list to the collection and adds their
        parameters to the parameterSet
        """
        if calc.GetName() in self:
            raise ValueError, "Calculation already has name %s" % str(calc.GetName())

        self[calc.GetName()] = calc 
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
                args = {'net': self[calc], 'vars': varsByCalc[calc],
                        'params': self.params}
                pypar.send((command, args), worker)

            # The master does his share here
            calc = calcs_to_do.pop()
            # We use the finally statement because we want to ensure that we
            #  *always* wait for replies from the workers, even if the master
            #  encounters an exception in his evaluation.
            try:
                results[calc] = self[calc].calculate(varsByCalc[calc], 
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
            calcPOrder = self[calcName].GetParameters().keys()
            self[calcName].CalculateSensitivity(varsByCalc[calcName], self.params)
            calcSensVals[calcName] = self[calcName].GetSensitivityResult(varsByCalc[calcName])
            calcVals[calcName] = self[calcName].GetResult(varsByCalc[calcName])

        return calcVals, calcSensVals

    def GetParameters(self):
        """
        Return a deep copy of the collections parameter KeyedList.
        """
        return copy.deepcopy(self.params)
