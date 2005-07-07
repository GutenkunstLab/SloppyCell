import sets, copy

from KeyedList import KeyedList

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
                    varsByCalc[calc][depVar].union_update(data[calc]\
                                                          [depVar].keys())

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
        self.set_shared_scale_factors(shared_sf)

    def SetName(self, name):
        self.name = name

    def GetName(self):
        return self.name

    def set_data(self, data):
        self.data = copy.copy(data)

    def UpdateData(self, newData):
        self.data.update(newData)

    def GetData(self):
        return self.data

    def set_fixed_scale_factors(self, fixed_sf):
        self.fixedScaleFactors = fixed_sf

    def set_shared_scale_factors(self, shared_sf):
        self.shared_sf = shared_sf

    def get_shared_scale_factors(self):
        return self.shared_sf

    def GetFixedScaleFactors(self):
        return self.fixedScaleFactors

    SetData = set_data
    SetFixedScaleFactors = set_fixed_scale_factors


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
        if params is not None:
            self.params.update(params)
	
        for (calcName, vars) in varsByCalc.items():
            self[calcName].Calculate(varsByCalc[calcName], self.params)

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
		
        for (calcName, vars) in varsByCalc.items():
            calcPOrder = self[calcName].GetParameters().keys()
            self[calcName].CalculateSensitivity(varsByCalc[calcName], self.params)

    def GetParameters(self):
        """
        GetParameterNamesInOrder() -> list

        Returns a list of parameter names, in the order they would be expected
        if a list of values were passed in.
        """
        return self.params
