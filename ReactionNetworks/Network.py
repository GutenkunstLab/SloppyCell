# This makes integer like one would expect mathematically, e.g. 1/2 = .5
#  rather than 0. This is obviously safer when loading other folks' models.
#  It does, however, cost about 10% performance on PC12.
# This is expected to become standard behavior around python 3.0
from __future__ import division

import copy
import exceptions
import sets
import sys
import types
import os

import scipy

from SloppyCell import _TEMP_DIR, KeyedList
if not os.path.isdir(_TEMP_DIR): 
    os.mkdir(_TEMP_DIR)

# This is the symbolic differentiation library
import symbolic
# We load a dictionary of previously-taken derivatives for efficiency
symbolic.loadDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))

import Integration
import Parsing
import Reactions
from Components import *
from Trajectory import Trajectory

# Expose function definitions that SBML wants to see.
log, log10 = scipy.log, scipy.log10
exp = scipy.exp
cos, sin, tan = scipy.cos, scipy.sin, scipy.tan
acos, asin, atan = scipy.arccos, scipy.arcsin, scipy.arctan
cosh, sinh, tanh = scipy.cosh, scipy.sinh, scipy.tanh
arccosh, arcsinh, arctanh = scipy.arccosh, scipy.arcsinh, scipy.arctanh
exponentiale, pi = scipy.e, scipy.pi

# This is the absolute path to the directory where Networks.py is. I use it
#  to find the xml files needed for the TeX output.
RXNNETS_DIR = os.path.dirname(__file__)

# Optional since it's x86-only.
try:
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# Optional since only TeXForm depends on them.
try:
    import Ft
    import libsbml
    HAVE_FT_AND_LIBSBML = True
except ImportError:
    HAVE_FT_AND_LIBSBML = False

class Network:

    def __init__(self, id, name=''):
        self.id, self.name = id, name

        self.functionDefinitions = KeyedList()
        self.reactions = KeyedList()
        self.assignmentRules = KeyedList()
        self.rateRules = KeyedList()
        self.events = KeyedList()

        # Variables is primary storage for all compartments, species, and 
        #  parameters. We build up 'cross reference' lists for convenience.
        self.variables = KeyedList()

        # The following are all 'cross reference' lists which are intended
        #  only to hold references to objects in self.variables
        #  (makeCrossReferences will fill them in)
        self.assignedVars = KeyedList()
        self.constantVars = KeyedList()
        self.optimizableVars = KeyedList()
        self.dynamicVars = KeyedList()
        self.compartments = KeyedList()
        self.parameters = KeyedList()
        self.species = KeyedList()

        # These are also cross references. complexEvents for events whose
        #  firing depends on variable values, and timeTriggeredEvents for
        #  events that fire at a give time.
        self.complexEvents = KeyedList()
        self.timeTriggeredEvents = KeyedList()

        # Add in the function definitions I've seen so far in SBML
        self._addStandardFunctionDefinitions()

        # Generating the functions in compile() is pretty slow, so dirty is 
        #  checked to see which ones need to be remade. (i.e. get_ddv_dt and 
        #  get_d2dv_ddvdt don't change if we add a new event, but do if we add
        #  a new variable.)
        self.dirty = {'ddv_dt': True, 'events': True}

        # Here are all the functions we dynamically generate. To add a new one,
        #  called, for example, fooFunc, you need to define make_fooFunc, which
        #  generates the python code for the function in fooFunc_functionBody.
        self.dynamicFunctions = ['get_ddv_dt', 'get_d2dv_ddvdt', 
                                 'get_d2dv_dovdt', 'get_eventValues', 
                                 'get_eventDerivs']

        # Should we get sensitivities via finite differences? (Faster, but less
        #  accurate.)
        self.fdSensitivities = False

    #
    # Methods used to build up a network
    #
    def add_variable(self, var):
        self._checkIdUniqueness(var.id)
        self.variables.setByKey(var.id, var)
        self.makeCrossReferences()

    def add_parameter(self, id, value=0.0, name='', is_constant=True, 
                      typical_value=None, is_optimizable=True):
        parameter = Parameter(id, value, name, is_constant, typical_value,
                              is_optimizable)
        self.add_variable(parameter)

    def add_compartment(self, id, size=1.0, name='', is_constant=True, 
                        typical_value=None, is_optimizable=False):
        """ Adds a compartment to the network. """
        compartment = Compartment(id, size, name, is_constant, typical_value,
                                  is_optimizable)
        self.add_variable(compartment)

    def add_species(self, id, compartment, initial_conc=None, 
                    name='', typical_value=None,
                    is_boundary_condition=False, is_constant=False, 
                    is_optimizable=False):
        species = Species(id, compartment, initial_conc, name, typical_value,
                          is_boundary_condition, is_constant, is_optimizable)
        self.add_variable(species)
        if species.id in self.dynamicVars.keys():
            self.dirty['ddv_dt'] = True

    def addEvent(self, *args, **kwargs):
        event = apply(Event, args, kwargs)
        self._checkIdUniqueness(event.id)
        self.events.setByKey(event.id, event)
        self.dirty['events'] = True
        self.makeCrossReferences()

    def addFunctionDefinition(self, *args, **kwargs):
        function = apply(FunctionDefinition, args, kwargs)
        self._checkIdUniqueness(function.id)
        self.functionDefinitions.setByKey(function.id, function)


    def addReaction(self, id, *args, **kwargs):
        # Reactions can be added by (1) passing in a string representing
        #  kinetic law, or (2) passing in a class already specifying the 
        #  kinetic law.
        # XXX: I'm a little unhappy with this because option (2) breaks the
        #      pattern that the first argument is the id
        if type(id) == types.StringType:
            rxn = apply(Reactions.Reaction, (id,) + args, kwargs)
        else:
            rxn = apply(id, args, kwargs)

        self._checkIdUniqueness(rxn.id)
        self.reactions.setByKey(rxn.id, rxn)
        self.dirty['ddv_dt'] = True

    def addAssignmentRule(self, variable, math):
        self.assignmentRules.setByKey(variable, math)
        self.makeCrossReferences()
        self.dirty['ddv_dt'] = True

    def add_rate_rule(self, lhs, rhs):
        self.rateRules.setByKey(lhs, rhs)
        self.makeCrossReferences()
        self.dirty['ddv_dt'] = True

    def _checkIdUniqueness(self, id):
        if id in self.variables.keys()\
           or id in self.reactions.keys()\
           or id in self.functionDefinitions.keys()\
           or id in self.events.keys()\
           or id == self.id:
            raise ValueError, ('The id %s is already in use!' % id)


    #
    # Methods to become a 'SloppyCell.Model'
    #

    def Calculate(self, vars, params = None):
        # Add in the times required by all the variables, then convert back to 
        #  a sorted list.
        # Make sure we start from t = 0
        t = sets.Set([0])
        for var, times in vars.items():
            t.union_update(times)
        t = list(t)
        t.sort()

        self.compile()

        self.trajectory = self.integrate(t, params, addTimes = True)

    def GetName(self):
        return self.id

    def GetParameters(self):
        return KeyedList([(var.id, var.value) for var in
                          self.optimizableVars.values()])

    def GetResult(self, vars):
        result = {}
        times = self.trajectory.timepoints
        for id in vars:
            traj = self.trajectory.getVariableTrajectory(id)
            result[id] = dict(zip(times, traj))

        return result

    #
    # The actual integration
    #
    def integrate(self, times, params = None,
                  returnEvents = False, addTimes = True,
                  rtol = None):
        if params is not None:
            self.setOptimizables(params)

        # We can add extra times to try to catch events better
        if addTimes:
            times = sets.Set(times)
            toAdd = scipy.linspace(min(times), max(times), 200)
            times.union_update(toAdd)
            times = list(times)
            times.sort()

        # If you ask for time = 0, we'll assume you want dynamic variable values
        #  reset.
        if times[0] == 0:
            self.resetDynamicVariables()

        # te, ye, and ie are the times, dynamic variable values, and indices
        #  of fired events
        t, oa, te, ye, ie = Integration.Integrate(self, times, rtol = rtol)
        
        indexKeys = self.dynamicVars.keys() + self.assignedVars.keys()
        keyToColumn = KeyedList(zip(indexKeys, range(len(indexKeys))))

	trajectory = Trajectory(self, keyToColumn)
        trajectory.appendFromODEINT(t, oa)

        if returnEvents:
            return trajectory, te, ye, ie
        else:
            return trajectory

    def integrateSensitivity(self, times, params = None,
                             returnEvents = False, addTimes = True,
                             rtol=None):
        if params is not None:
            self.setOptimizables(params)

        # We can add extra times to try to catch events better
        if addTimes:
            times = sets.Set(times)
            toAdd = scipy.linspace(min(times), max(times), 200)
            times.union_update(toAdd)
            times = list(times)
            times.sort()

        # If you ask for time = 0, we'll assume you want dynamic variable values
        #  reset.
        if times[0] == 0:
            self.resetDynamicVariables()

        nOv, nDv, nAv = len(self.optimizableVars), len(self.dynamicVars),\
                len(self.assignedVars)

        if not self.fdSensitivities:
            t, oa, te, ye, ie = Integration.Integrate_Ddv_Dov(self, times,
                                                              rtol=rtol)
	else:
            # We do it by finite differencing
            oa = scipy.zeros((len(times), nDv*(nOv+1)), scipy.Float)
            t, oaInitial, te, ye, ie = Integration.Integrate(self, times,
                                                             rtol=rtol)
            oa[:,0:nDv] = oaInitial
            for i, (id, var) in enumerate(self.optimizableVars.items()) :
                saved = var.value

                stepsize = 1.0e-6*abs(saved) + 1.0e-12
                self.setInitialVariableValue(id, saved + stepsize)

	        self.resetDynamicVariables()
                tStepped, oaStepped, teStepped, yeStepped, ieStepped = \
                        Integration.Integrate(self, times, rtol = rtol)
                #XXX: Should pull values corresponding to times out of tStepped,
                # oaStepped since param change might mean events fire at 
                # different times or not at all.
		deriv = (oaStepped - oaInitial)/stepsize
		self.setInitialVariableValue(id, saved)

                oa[:, nDv + i*nDv : nDv + (i+1)*nDv] = deriv

	alldvs = self.dynamicVars.keys() + self.assignedVars.keys()
        keyToColumnSensNames = alldvs + [(cname, pname) for
                                         pname in self.optimizableVars.keys()
                                         for cname in alldvs]
        keyToColumnSens = [(name, index) for index, name in
                           enumerate(keyToColumnSensNames)]
	keyToColumnSens = KeyedList(keyToColumnSens)

	ddv_dpTrajectory = Trajectory(self, keyToColumnSens)
	ddv_dpTrajectory.appendSensFromODEINT(t, oa)

        if returnEvents:
            return ddv_dpTrajectory, te, ye, ie
        else:
            return ddv_dpTrajectory

    def FindFixedPoint(self, dv0 = None):
        if dv0 is None:
            dv0 = self.getDynamicVarValues()

        ddv_dtFromLogs = lambda logDV: self.get_ddv_dt(scipy.exp(logDV), 0)
        fprime = lambda logDV: self.get_d2dv_ddvdt(scipy.exp(logDV), 0)\
                *scipy.exp(logDV)

        dvFixed = scipy.optimize.fsolve(ddv_dtFromLogs, x0 = scipy.log(dv0),
                                        fprime = fprime, col_deriv = True)

        return scipy.exp(dvFixed)

    def PrintFakeData(self, times, dataIds = None,
                      fixedScaleFactors = False, fileName = None,
                      sigmaLambda = lambda y: 0*y + 1):
        if dataIds == None:
            dataIds = self.dynamicVars.keys() + self.assignedVars.keys()
        if fileName == None:
            fileName = "fakedata_%s.py" % self.id
        
        traj, te, ye, ie = self.integrate(times, addTimes = False, returnEvents = True)

        output = file(fileName, 'w')
        output.write("import SloppyCell.Collections as Collections\n")
        output.write("fakedata = Collections.Experiment('fakeExpt')\n")
        output.write("fakedata.longname = 'Data generated directly from network %s'\n" % self.id)

        if fixedScaleFactors:
            output.write('fakedata.SetFixedScaleFactors({%s})\n' %
                         ',\n'.join(["'%s': 1.0" % id for id in dataIds])
                         )

	output.write("fakedata.SetData({'%s':{\n" % self.id)
        for dataId in dataIds:
            output.write("\t\t\t'%s':{\n" % dataId)
            y = traj.getVariableTrajectory(dataId)
            sigma = sigmaLambda(y)
            for time, value, sig in zip(traj.timepoints, y, sigma):
                if time not in te:
                    output.write("\t\t\t\t%f: (%f, %f),\n" % (time, value, sig))
            output.write("\t\t\t},##end %s data\n" % dataId)
            
        output.write("\t\t},##end calc data\n")
        output.write("\t})##end fakedata.SetData\n")
        output.write("fakeDataColl = Collections.ExperimentCollection([fakedata])\n")
        output.close()


    #
    # Methods to get and set object properties
    #

    def setInitialVariableValue(self, id, value):
        var = self.variables.getByKey(id)
        var.initialValue = value
        if var.isConstant:
            var.value = value
        self.constantVarValues = [var.value for var in
                                  self.constantVars.values()]

    def setTypicalVariableValue(self, id, value):
        self.variables.getByKey(id).typicalValue = value

    def setOptimizables(self, params):
        # Sets all the optimizable variables from a passed-in KeyedList
        for id, value in params.items():
            self.setInitialVariableValue(id, value)

    def getInitialVariableValue(self, id):
        return self.variables.getByKey(id).initialValue

    def getDynamicVarValues(self):
        return [var.value for var in self.dynamicVars.values()]

    def resetDynamicVariables(self):
        # Resets all dynamical variables to their initial values. This is
        #  a little complex because the initial value may be a function of
        #  other values. Thus we skip those ones on the first pass through.
        pending = []
        for id, var in self.dynamicVars.items():
            if type(var.initialValue) == types.StringType:
                pending.append((id, var))
            else:
                var.value = var.initialValue

        self.updateAssignedVars(time = 0)

        for id, var in pending:
            rhs = self.substituteFunctionDefinitions(var.initialValue)
            rhs = self.substituteVariableNames(rhs)
            self.dynamicVars.getByKey(id).value = eval(rhs)

    def updateVariablesFromDynamicVars(self, values, time):
        for ii in range(len(self.dynamicVars)):
            self.dynamicVars[ii].value = values[ii]
        self.updateAssignedVars(time)

    def updateAssignedVars(self, time):
        for id, rhs in self.assignmentRules.items():
            rhs = self.substituteFunctionDefinitions(rhs)
            rhs = self.substituteVariableNames(rhs)
            self.assignedVars.getByKey(id).value = eval(rhs)

    #
    # Generate the differential equations and functions to calculate them.
    #

    def makeDiffEqRHS(self):
        # Start with all right-hand-sides set to '0'
        self.diffEqRHS = KeyedList(zip(self.dynamicVars.keys(), 
                                       ['0'] * len(self.dynamicVars)))

        for id, rxn in self.reactions.items():
            rateExpr = rxn.kineticLaw
            for reactantId, dReactant in rxn.stoichiometry.items():
                if self.variables.getByKey(reactantId).is_boundary_condition:
                    # Variables that are boundary conditions aren't modified by
                    #  reactions, so we move on to the next reactant.
                    continue

                reactantIndex = self.dynamicVars.indexByKey(reactantId)
                if (type(dReactant) == types.IntType
                    or type(dReactant) == types.FloatType)\
                   and dReactant != 0:
                    if dReactant == 1:
                        self.diffEqRHS[reactantIndex] += ' + ( %s )' % rateExpr
                    elif dReactant == -1:
                        self.diffEqRHS[reactantIndex] += ' - ( %s )' % rateExpr
                    else:
                        self.diffEqRHS[reactantIndex] += ' + %f * ( %s )'\
                                % (dReactant, rateExpr)

                elif type(dReactant) == types.StringType:
                    self.diffEqRHS[reactantIndex] += ' + ( %s ) * ( %s )'\
                            % (dReactant, rateExpr)

        # Strip off the initial '0 +/- '
        for id, rhs in self.diffEqRHS.items():
            if len(rhs) > 4:
                if rhs[:4] == '0 + ':
                    self.diffEqRHS.setByKey(id, rhs[4:])
                elif rhs[:4] == '0 - ':
                    self.diffEqRHS.setByKey(id, rhs[2:])

        # Handle the rate rules
        for id, rule in self.rateRules.items():
            self.diffEqRHS.setByKey(id, rule)

    def make_get_ddv_dt(self):
        self.ddv_dt = scipy.zeros(len(self.dynamicVars), scipy.Float)

        functionBody = 'def get_ddv_dt(self, dynamicVars, time):\n\t'
        functionBody += 'ddv_dt = self.ddv_dt\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)
        functionBody += '\n\t'

        for ii, (id, var) in enumerate(self.dynamicVars.items()):
            rhs = self.diffEqRHS.getByKey(id)
            if rhs != '0':
                rhs = self.substituteFunctionDefinitions(rhs)
                rhs = self.substituteConstantVariableValues(rhs)
                functionBody += '# Total derivative of %s wrt time\n\t' % (id)
                functionBody += 'ddv_dt[%i] = %s\n\t' % (ii, rhs)

        functionBody += '\n\n\treturn ddv_dt'

        f = file(os.path.join(_TEMP_DIR, 'get_ddv_dt.py'), 'w')
        print >> f, functionBody
        f.close()
        self.get_ddv_dt_functionBody = functionBody
        self.get_ddv_dt = None
        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))

    def make_get_d2dv_ddvdt(self):
        self.d2dv_ddvdt = scipy.zeros((len(self.dynamicVars),
                                     len(self.dynamicVars)), scipy.Float)

        functionBody = 'def get_d2dv_ddvdt(self, dynamicVars, time):\n\t'
        functionBody += 'd2dv_ddvdt = self.d2dv_ddvdt\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)
        functionBody += '\n\t'

        for rhsIndex, rhsId in enumerate(self.dynamicVars.keys()):
            rhs = self.diffEqRHS.getByKey(rhsId)
            rhs = self.substituteFunctionDefinitions(rhs)
            rhs = self.substituteConstantVariableValues(rhs)
            # Take all derivatives of the rhs with respect to other 
            #  dynamic variables
            for wrtIndex, wrtId in enumerate(self.dynamicVars.keys()):
                deriv = self.takeDerivative(rhs, wrtId)
                if deriv != '0':
                    functionBody += '# Partial derivative of Ddv_Dt for %s wrt %s\n\t' % (rhsId, wrtId)
                    functionBody += 'd2dv_ddvdt[%i, %i] = %s\n\t' % \
                            (wrtIndex, rhsIndex, deriv)

        functionBody += '\n\treturn d2dv_ddvdt'

        f = file(os.path.join(_TEMP_DIR, 'get_d2dv_ddvdt.py'), 'w')
        print >> f, functionBody
        f.close()
        self.get_d2dv_ddvdt_functionBody = functionBody
        self.get_d2dv_ddvdt = None
        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))

    def make_get_d2dv_dovdt(self):
        self.d2dv_dovdt = scipy.zeros((len(self.dynamicVars),
                                       len(self.optimizableVars)), scipy.Float)

        functionBody = 'def get_d2dv_dovdt(self, dynamicVars, time, indices = None):\n\t'
        functionBody += 'd2dv_dovdt = self.d2dv_dovdt\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)
        functionBody += '\n\t'

        for wrtIndex, wrtId in enumerate(self.optimizableVars.keys()):
            functionBody += 'if indices is None or %i in indices:\n\t\t' % wrtIndex
            derivWritten = False
            for rhsIndex, rhsId in enumerate(self.dynamicVars.keys()):
                rhs = self.diffEqRHS.getByKey(rhsId)
                rhs = self.substituteFunctionDefinitions(rhs)
                rhs = self.substituteConstantVariableValues(rhs)
                deriv = self.takeDerivative(rhs, wrtId)
                if deriv != '0':
                    functionBody += '# Partial derivative of Ddv_Dt for %s wrt %s\n\t\t' % (rhsId, wrtId)
                    functionBody += 'd2dv_dovdt[%i, %i] = %s\n\t\t' % \
                            (rhsIndex, wrtIndex, deriv)
                    derivWritten = True
            if derivWritten == False:
                functionBody += 'pass\n\t\t'
            functionBody += '\n\t'

        functionBody += '\n\treturn d2dv_dovdt'

        f = file(os.path.join(_TEMP_DIR, 'get_d2dv_dovdt.py'), 'w')
        print >> f, functionBody
        f.close()
        self.get_d2dv_dovdt_functionBody = functionBody
        self.get_d2dv_dvcdt = None
        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))


    #
    # Methods involved in events
    #

    def fireEvent(self, event, dynamicVarValues, time):
        self.updateVariablesFromDynamicVars(dynamicVarValues, time)

        delay = event.delay
        # If the delay is a math expression, calculate its value
        if type(delay) == types.StringType:
            delay = self.substituteFunctionDefinitions(delay)
            delay = self.substituteVariableNames(delay)
            delay = eval(delay)

        executionTime = round(time + delay, 8)
        event.newValues[executionTime] = {}
        # Calculate the values that will get assigned to the variables
        for id, rhs in event.eventAssignments.items():
            if type(rhs) == types.StringType:
                rhs = self.substituteFunctionDefinitions(rhs)
                erhs = self.substituteVariableNames(rhs)
                event.newValues[executionTime][id] = eval(erhs)
            else:
                event.newValues[executionTime][id] = rhs

        return delay

    def fireEventAndUpdateDdv_Dov(self, event, values, time):
        # XXX: This a big mess, but it works. Should refactor.
        self.updateVariablesFromDynamicVars(values, time)

        nDV, nOV = len(self.dynamicVars), len(self.optimizableVars)
        delay, trigger = event.delay, event.trigger

        # If the delay is a math expression, calculate its value
        if type(delay) == types.StringType:
            if len(Parsing.extractVariablesFromString(delay)) > 0:
                raise exceptions.NotImplementedError, "We don't support math form delays in sensitivity! (Yet)"
            else:
                delay = eval(delay)

        executionTime = round(time + delay, 8)
        event.newValues[executionTime] = {}

        ddv_dt = self.get_ddv_dt(values[:nDV], time)
        newValues = copy.copy(values[:nDV])
        for lhsId, rhs in event.eventAssignments.items():
            lhsIndex = self.dynamicVars.indexByKey(lhsId)

            # Everything below should work if the rhs is a constant, as long
            #  as we convert it to a string first.
            rhs = str(rhs)
            rhs = self.substituteFunctionDefinitions(rhs)

            # Compute and store the new value of the lhs
            newVal = eval(self.substituteVariableNames(rhs))
            event.newValues[executionTime][lhsId] = newVal
            newValues[lhsIndex] = newVal

        trigger = self.substituteFunctionDefinitions(event.trigger)
        # Compute the derivatives of our firing times wrt all our variables
        # dTf_dov = -(sum_dv dTrigger_ddv ddv_dov - dTrigger_dov)/(sum_dv dTrigger_ddv ddv_dt)
        # XXX: Handle trigger has explicit time dependence.
        if event.timeTriggered:
            dTf_dov = dict(zip(self.optimizableVars.keys(), [0] * nOV))
        else:
            if 'time' in Parsing.extractVariablesFromString(trigger):
                raise exceptions.NotImplementedError, "We don't support explicit time dependence in complex event triggers for sensitivity! (Yet) Tigger was: %s" % trigger
            dtrigger_ddvValue = {}
            for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                dtrigger_ddv = self.takeDerivative(trigger, dynId)
                dtrigger_ddv = self.substituteVariableNames(dtrigger_ddv)
                dtrigger_ddvValue[dynId] = eval(dtrigger_ddv)

            dtrigger_dovValue = {}
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                dtrigger_dov = self.takeDerivative(trigger, optId)
                dtrigger_dov = self.substituteVariableNames(dtrigger_dov)
                dtrigger_dovValue[optId] = eval(dtrigger_dov)

            dTf_dov = {}
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                numerator, denominator = 0, 0
                for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                    indexInOA = dvIndex + (ovIndex + 1)*nDV
                    denominator += dtrigger_ddvValue[dynId] *\
                            self.get_ddv_dt(values[:nDV], time)[dvIndex]
                    numerator += dtrigger_ddvValue[dynId] * values[indexInOA]

                numerator += dtrigger_dovValue[optId]
                dTf_dov[optId] = -numerator/denominator

        for lhsIndex, lhsId in enumerate(self.dynamicVars.keys()):
            rhs = event.eventAssignments.get(lhsId, lhsId)
            # Everything below should work if the rhs is a constant, as long
            #  as we convert it to a string first.
            rhs = str(rhs)
            rhs = self.substituteFunctionDefinitions(rhs)

            drhs_ddvValue = {}
            for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                drhs_ddv = self.takeDerivative(rhs, dynId)
                drhs_ddv = self.substituteVariableNames(drhs_ddv)
                drhs_ddvValue[dynId] = eval(drhs_ddv)

            drhs_dovValue = {}
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                drhs_dov = self.takeDerivative(rhs, optId)
                drhs_dov = self.substituteVariableNames(drhs_dov)
                drhs_dovValue[optId] = eval(drhs_dov)

            drhs_dtime = self.takeDerivative(rhs, 'time')
            drhs_dtime = self.substituteVariableNames(drhs_dtime)
            drhs_dtimeValue = eval(drhs_dtime)

            # Calculate the perturbations to the sensitivities
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                    event.newValues[executionTime].setdefault((lhsId, optId), 0)
                    indexInOA = dvIndex + (ovIndex + 1)*nDV
                    #dv_ddv * ddv_dov
                    event.newValues[executionTime][(lhsId, optId)] += \
                            drhs_ddvValue[dynId] * values[indexInOA]
                    #dv_ddv * ddv_dt * dTf_dov
                    event.newValues[executionTime][(lhsId, optId)] += \
                            self.get_ddv_dt(values[:nDV], time)[dvIndex]\
                            * drhs_ddvValue[dynId] * dTf_dov[optId]

                event.newValues[executionTime].setdefault((lhsId, optId), 0)

                #dv_dov
                event.newValues[executionTime][(lhsId, optId)] += drhs_dovValue[optId]

                #dv_dtime
                event.newValues[executionTime][(lhsId, optId)] += \
                        drhs_dtimeValue * dTf_dov[optId]

                # dv_dt * dTf_dov
                # note that dc_dt * dTf_dov cancels out
                event.newValues[executionTime][(lhsId, optId)] += \
                        -self.get_ddv_dt(newValues, time)[lhsIndex] * dTf_dov[optId]

        return delay

    def executeEvent(self, event, dynamicVarValues, time):
        for id, value in event.newValues[round(time, 8)].items():
            if type(id) == types.StringType:
                dynamicVarValues[self.dynamicVars.indexByKey(id)] = value

        return dynamicVarValues

    def executeEventAndUpdateDdv_Dov(self, event, inValues, time):
        values = copy.copy(inValues)
        # Copy out the new values for the dynamic variables
        for id, value in event.newValues[round(time, 8)].items():
            if type(id) == types.TupleType:
                nDV = len(self.dynamicVars)
                dvIndex = self.dynamicVars.indexByKey(id[0])
                ovIndex = self.optimizableVars.indexByKey(id[1])
                values[dvIndex + (ovIndex + 1)*nDV] = value

            else:
                values[self.dynamicVars.indexByKey(id)] = value

        return values


    def make_get_eventValues(self):
        self._eventValues = scipy.zeros(len(self.complexEvents), scipy.Float)
        self._eventTerminals = scipy.zeros(len(self.complexEvents))
        self._eventDirections = scipy.zeros(len(self.complexEvents))

        functionBody = 'def get_eventValues(self, dynamicVars, time):\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)

        for ii, event in enumerate(self.complexEvents.values()):
            rhs = self.substituteFunctionDefinitions(event.trigger)
            functionBody += 'self._eventValues[%i] = %s\n\t' % (ii, rhs)
            functionBody += 'self._eventTerminals[%i] = %s\n\t' % \
                    (ii, event.isTerminal)
            functionBody += 'self._eventDirections[%i] = 1\n\t' % ii

        functionBody += '\n\treturn self._eventValues, self._eventTerminals, self._eventDirections'

        f = file(os.path.join(_TEMP_DIR, 'get_eventValues.py'), 'w')
        print >> f, functionBody
        f.close()
        self.get_eventValues_functionBody = functionBody
        self.get_eventValues = None

    def make_get_eventDerivs(self):
        self._eventDerivValues = scipy.zeros(len(self.complexEvents),
                                             scipy.Float)
        
        functionBody = 'def get_eventDerivs(self, dynamicVars, ddv_dt, time):\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)

        for ii, event in enumerate(self.complexEvents.values()):
            trigger = self.substituteFunctionDefinitions(event.trigger)
            rhs = ''
            for id in self.diffEqRHS.keys():
                # We use the chain rule to get derivatives wrt time.
                deriv = self.takeDerivative(trigger, id)
                if deriv != '0':
                    rhs += ' + (%s) * ddv_dt[%i]' % \
                            (deriv, self.dynamicVars.indexByKey(id))

            # We need to include the partial derivative wrt time.
            deriv = self.takeDerivative(trigger, 'time')
            if deriv != '0':
                rhs += ' + (%s)' % deriv

            rhs = rhs[3:]
            functionBody += 'self._eventDerivValues[%i] = %s\n\t'\
                    % (ii, rhs)

        functionBody += '\n\treturn self._eventDerivValues'

        f = file(os.path.join(_TEMP_DIR, 'get_eventDerivs.py'), 'w')
        print >> f, functionBody
        f.close()
        self.get_eventDerivs_functionBody = functionBody
        self.get_eventDerivs = None

    #
    # Internally useful things
    #

    def makeCrossReferences(self):
        self.assignedVars = KeyedList()
        self.constantVars = KeyedList()
        self.optimizableVars = KeyedList()
        self.dynamicVars = KeyedList()

        self.compartments = KeyedList()
        self.parameters = KeyedList()
        self.species = KeyedList()
        mapping = {Compartment: self.compartments,
                   Parameter: self.parameters,
                   Species: self.species}

        for id, var in self.variables.items():
            mapping[var.__class__].setByKey(id, var)
            if var.isConstant:
                self.constantVars.setByKey(id, var)
                if var.isOptimizable:
                    self.optimizableVars.setByKey(id, var)
            elif id in self.assignmentRules.keys():
                self.assignedVars.setByKey(id, var)
            else:
                self.dynamicVars.setByKey(id, var)

        self.constantVarValues = [var.value for var in 
                                  self.constantVars.values()]

        self.complexEvents = KeyedList()
        self.timeTriggeredEvents = KeyedList()
        for id, event in self.events.items():
            if event.timeTriggered:
                self.timeTriggeredEvents.setByKey(id, event)
            else:
                self.complexEvents.setByKey(id, event)

    def compile(self):
        # XXX: Should really figure out a better way to do this.
        if self.dirty['ddv_dt']:
            self.makeDiffEqRHS()
            self.make_get_ddv_dt()
            self.make_get_d2dv_ddvdt()
            self.make_get_d2dv_dovdt()
            self.dirty['ddv_dt'] = False

        if self.dirty['events']:
            self.make_get_eventValues()
            self.make_get_eventDerivs()
            self.dirty['events'] = False

        # Create all the functions listed in self.dynamicFunctions, if they
        #  don't already exist.
        for func in self.dynamicFunctions:
            if getattr(self, func, None) is None:
                exec getattr(self, func+'_functionBody')
                setattr(self, func, 
                        types.MethodType(locals()[func], self, Network))
                if HAVE_PSYCO:
                    psyco.bind(getattr(self, func))


    def addAssignmentRulesToFunctionBody(self, functionBody):
        functionBody += 'constantVarValues = self.constantVarValues\n\n\t'

        for ii, id in enumerate(self.constantVars.keys()):
            functionBody += '%s = constantVarValues[%i]\n\t' % (id, ii)

        functionBody += '\n\t'
        for ii, id in enumerate(self.dynamicVars.keys()):
            functionBody += '%s = dynamicVars[%i]\n\t' % (id, ii)

        for variable, math in self.assignmentRules.items():
            rhs = self.substituteFunctionDefinitions(math)
            functionBody += variable + ' = ' + rhs + '\n\t'

        return functionBody

    def takeDerivative(self, input, wrt):
        """
        Take the derivative of a math expression wrt a given variable id.

        Does the chain rule through assigned variables.
        """
        output = symbolic.Diff(input, wrt)

        # What other assigned variables does input depend on?
        otherVars = Parsing.extractVariablesFromString(input)
        otherVars.difference_update([wrt])
        otherVars.intersection_update(self.assignedVars.keys())
        # Do the chain rule for those variables
        if id in self.assignedVars.keys():
            rule = self.assignmentRules.getByKey(id)
            rule = self.substituteFunctionDefinitions(rule)
            d2 = self.takeDerivative(rule, wrt)
            if d2 != '0':
                d = symbolic.Diff(input, id)
                output += ' + %s *(%s)' % (d, d2)

        return output

    def substituteConstantVariableValues(self, input):
        for id, var in self.constantVars.items():
            if id not in self.optimizableVars.keys():
                input = Parsing.substituteVariableNamesInString(input, id,
                                                                str(var.value)) 
        return input

    def substituteVariableNames(self, input):
        for id in Parsing.extractVariablesFromString(input):
            if id != 'time':
                mapping = "self.variables.getByKey('%s').value" % id
                input = Parsing.substituteVariableNamesInString(input, id,
                                                                mapping) 
        return input

    def substituteFunctionDefinitions(self, input):
        reversedFunctions = self.functionDefinitions.items()
        reversedFunctions.reverse()
        functionsUsed = Parsing.extractFunctionsFromString(input)
        for id, function in reversedFunctions:
            if id in functionsUsed:
                input = Parsing.substituteFunctionIntoString(input, id, 
                                                             function.variables,
                                                             function.math)
                functionsUsed = Parsing.extractFunctionsFromString(input)

        standardIds = sets.Set(self._standardFuncDefs.keys())
        standardUsed = standardIds.intersection(functionsUsed)
        while len(standardUsed) > 0:
            for id in standardUsed:
                function = self._standardFuncDefs.getByKey(id)
                input = Parsing.substituteFunctionIntoString(input, id, 
                                                             function.variables,
                                                             function.math)

            functionsUsed = Parsing.extractFunctionsFromString(input)
            standardUsed = standardIds.intersection(functionsUsed)


        return input

    def TeXForm(self):
        """
        Return a dict of TeX-formatted network eqns, indexed by variable id.

        The flow here is a little complicated. We use libsbml to go from the
        equations in python syntax to MathML. Then we use 4Suite and a set of
        XML style sheets to convert the MathML to TeX.
        """
        if not HAVE_FT_AND_LIBSBML:
            raise exceptions.Warning, 'TeXForm not available because 4Suite library or libsbml is not installed.'
            return

        # I don't fully understand the usage of Ft, so much is stolen from
        # http://uche.ogbuji.net/tech/akara/nodes/2003-01-01/python-xslt
        from Ft.Xml.Xslt import Processor
        processor = Processor.Processor()
        from Ft.Xml import InputSource

        # This is the absolute path to the stylesheets we use
        transformFileName = os.path.join(RXNNETS_DIR, 'xsltml_2.0', 
                                         'mmltex.xsl')
        transformFile = file(transformFileName, 'r')

        # These substitutions munge Windows paths into appropriate URIs
        transformURI = transformFileName
        transformURI = transformURI.replace(':', '|')
        transformURI = transformURI.replace('\\', '/')
        transformURI = 'file://' + transformFileName

        transform = InputSource.DefaultFactory.fromStream(transformFile, 
                                                          transformURI)
        processor.appendStylesheet(transform)

        # libsbml wants the MathML embedded in a MathMLDocument
        mmlDoc = libsbml.MathMLDocument()
        eqns = KeyedList()
        for id, rhs in self.diffEqRHS.items() + self.assignmentRules.items():
            # Convert powers from python notation to what libsbml expects
            rhs = rhs.replace('**', '^')

            mathAST = libsbml.parseFormula(rhs)
            mmlDoc.setMath(mathAST)
            sourceString = libsbml.writeMathMLToString(mmlDoc)
            source = InputSource.DefaultFactory.fromString(sourceString, 
                                                           'file://bogus')
            # Finally, we get the TeX-formatted RHS
            texRHS = processor.run(source)

            # Add in the LHS of eqn
            if id in self.assignmentRules.keys():
                line = '\\mathrm{%s} &=& %s' % (id, texRHS[1:-1])
            else:
                line = '{d\\mathrm{%s} \\over dt} &=& %s' % (id, texRHS[1:-1])

            # Put \times between variables for clarity
            line = line.replace('}\\mathrm', '}\\times\\mathrm')

            # Replace variable ids with their names
            for id2 in self.variables.keys():
                # Wrap species in brackets
                if id2 in self.species.keys():
                    line = line.replace('\\mathrm{%s}' % id2,
                                        '\\left[\\mathrm{%s}\\right]' % id2)

                name = self.variables.getByKey(id2).name
                if name:
                    line = line.replace('\\mathrm{%s}' % id2,
                                        '\\mathrm{%s}' % name)
                else:
                    # If a parameter doesn't have a name, kludge together
                    #  a TeX name by subscripting anything after the first _
                    if id2 in self.parameters.keys() and id2.count('_') > 0:
                        sp = id2.split('_')
                        newId = '%s_{%s}' % (sp[0], ''.join(sp[1:]))
                        line = line.replace('\\mathrm{%s}' % id2,
                                            '\\mathrm{%s}' % newId)

            eqns.setByKey(id, line + '\\\\')

        return eqns

    def TeXOutputToFile(self, filename):
        """Output TeX-formatted network eqns to a file."""
        eqns = self.TeXForm()
        if eqns is None: # Can happen if we don't have Ft or libsbml
            return

	output = open(filename, 'w')
	output.write('\\documentclass[12pt]{article}\n')
        output.write('\\usepackage[pdftex, landscape]{geometry}\n')
	output.write('\\begin{document}\n')
	output.write('\\begin{eqnarray*}\n')
        for equation in eqns.values():
            output.write(equation)
            output.write('\n')
	output.write('\\end{eqnarray*}\n')
	output.write('\\end{document}\n')
	output.close()

    def __getstate__(self):
        odict = copy.copy(self.__dict__)

        for funcName in self.dynamicFunctions:
            odict[funcName] = None

        return odict

    def __setstate__(self, newdict):
        self.__dict__.update(newdict)
        self.makeCrossReferences()

    def _addStandardFunctionDefinitions(self):
        standardFuncDefs = {
                            'gt':(['x', 'y'],  'x - y'),
                            'geq':(['x', 'y'],  'x - y'),
                            'lt':(['x', 'y'],  'y - x'),
                            'leq':(['x', 'y'],  'y - x'),
                            'pow':(['x', 'n'],  'x**n'),
                            'root':(['n', 'x'],  'x**(1./n)'),
                            'sqrt':(['x'],  'x**(.5)'),
                            'cot':(['x'],  '1./tan(x)'),
                            'arccot':(['x'],  'atan(1./x)'),
                            'coth':(['x'],  '1./tanh(x)'),
                            'arccoth':(['x'],  'arctanh(1./x)'),
                            'csc':(['x'],  '1./sin(x)'),
                            'arccsc':(['x'],  'asin(1./x)'),
                            'csch':(['x'],  '1./sinh(x)'),
                            'arccsch':(['x'],  'arcsinh(1./x)'),
                            'sec':(['x'],  '1./cos(x)'),
                            'arcsec':(['x'],  'acos(1./x)'),
                            'sech':(['x'],  '1./cosh(x)'),
                            'arcsech':(['x'],  'arccosh(1./x)'),
                            }

        self._standardFuncDefs = KeyedList()
        for id, (vars, math) in standardFuncDefs.items():
            function = FunctionDefinition(id, vars, math)
            self._checkIdUniqueness(function.id)
            self._standardFuncDefs.setByKey(function.id, function)

        # XXX: SBML takes log(b, x) to mean log base b of x. We should
        #      modify Parsing.substituteFunctionIntoString to differentiate
        #      between functions with different numbers of arguments
        #self.addFunctionDefinition('log', ['b, x'],  'log(x)/log(b)')

    def addCompartment(self, id, size=1.0, name='', isConstant=True, 
                       typicalValue=False, isOptimizable=False):
        compartment = Compartment(id, size, name, isConstant, typicalValue,
                                  isOptimizable)
        self.addVariable(compartment)
        #raise exceptions.DeprecationWarning('Method addCompartment is deprecated, use add_compartment instead.')

    addVariable = add_variable
    
    def addSpecies(self, id, compartment, initialConcentration=None,
                   name='', typicalValue=None, is_boundary_condition=False,
                   isConstant=False, isOptimizable=False):
        self.add_species(id, compartment, initialConcentration, name,
                         typicalValue, is_boundary_condition, isConstant,
                         isOptimizable)

    def addParameter(self, id, value = 0.0, 
                     typicalValue = None, name = '',
                     isConstant = True, isOptimizable = True):
        self.add_parameter(id, value, name, isConstant, typicalValue,
                           isOptimizable)

    addRateRule = add_rate_rule
