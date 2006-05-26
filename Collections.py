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

    def GetResults(self, requestedByCalc):
        """
        GetResultsByCalc(requestedByCalc, params) -> dictionary

        Given requestedByCalc, a dictionary of the form:
         dict[calc name][dep var][ind var], and a set of parameters
        (dictionary or appropriately ordered list), returns a dictionary of
        results. The dictionary is of the form:
         dictionary[calculation][dependent variables][independent variable]
          -> result
        """
        calcVals = {}
        for (calcName, requested) in requestedByCalc.items():
            calcVals[calcName] = self[calcName].GetResult(requested)

        return calcVals

    def Calculate(self, varsByCalc, params = None):
        # This function has been parallelized. Each worker (including the root
        # node) does a subset of the calculations.
        my_rank, num_procs = SloppyCell.my_rank, SloppyCell.num_procs

        if params is not None:
            self.params.update(params)

        results = {}
        # Divide up the jobs amongst the various workers
        calcs_for_me = varsByCalc.keys()[my_rank::num_procs]

        # We need to be quite careful here about exceptions.
        # If a worker encounters a SloppyCellException, we assume it's probably
        #  an integration error. So we send it to the master and let it deal
        #  with it how it likes.
        # If a worker encounters any other exception, it's probably a sign of
        #  a bug in the code, so we don't catch it here and instead let 
        #  Model.MasterSwitch deal with it.
        for calc in calcs_for_me:
            if my_rank != 0:
                try:
                    self[calc].Calculate(varsByCalc[calc], self.params)
                except Utility.SloppyCellException, X:
                    logger.debug('SloppyCellException caught by processor %i'
                                 % my_rank)
                    results[calc] = X
                else:
                    results[calc] = self[calc].GetResult(varsByCalc[calc])
            else:
                # We don't want to catch and re-raise exceptions on the master
                #  node because it makes tracebacks much less useful.
                self[calc].Calculate(varsByCalc[calc], self.params)
                results[calc] = self[calc].GetResult(varsByCalc[calc])

        # If we're running parallel. Collect our results on the root node.
        if num_procs > 1 and my_rank == 0:
            for worker in range(1, num_procs):
                worker_results = pypar.receive(worker)
                results.update(worker_results)
            # The root node checks if the workers had any SloppyCellExceptions 
            #  and reraises them. It is very important that this be done *after*
            #  receiving results from every worker.
            for result in results.values():
                if isinstance(result, Utility.SloppyCellException):
                    raise result
        elif my_rank > 0:
            pypar.send(results, 0)

        return results 

    def GetSensitivityResults(self, requestedByCalc):
        """
        Given requestedByCalc, a dictionary of the form:
         dict[calc name][dep var][ind var], and a set of parameters
        (dictionary or appropriately ordered list), returns a dictionary of
        results. The dictionary is of the form:
         dictionary[calculation][dependent variables][independent variable]
	 [parameter]
          -> result
        """
        calcSensVals = {}
        for (calcName, requested) in requestedByCalc.items():
            calcSensVals[calcName] = self[calcName].GetSensitivityResult(requested)

        return calcSensVals


    def CalculateSensitivity(self, varsByCalc, params = None):
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
