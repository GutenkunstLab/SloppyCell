# This makes integer like one would expect mathematically, e.g. 1/2 = .5
#  rather than 0. This is obviously safer when loading other folks' models.
#  It does, however, cost about 10% performance on PC12.
# This is expected to become standard behavior around python 3.0
from __future__ import division

import copy
import exceptions
import sets
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
import IO
import Parsing
import Reactions
import SloppyCell.Collections as Collections
import matplotlib.mlab as mlab

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

# Optional since it's x86-only.
try:
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

try:
    import Dynamics
    HAVE_DYNAMICS = True
except ImportError:
    HAVE_DYNAMICS = False

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
        #  (_makeCrossReferences will fill them in)
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

        # Here are all the functions we dynamically generate. To add a new one,
        #  called, for example, fooFunc, you need to define _make_fooFunc, which
        #  generates the python code for the function in fooFunc_functionBody.
        self._dynamic_structure_funcs = ['get_ddv_dt', 'get_d2dv_ddvdt', 
                                         'get_d2dv_dovdt']
        self._dynamic_event_funcs = ['get_eventValues', 'get_eventDerivs', 
                                     'root_func']

        # Should we get sensitivities via finite differences? (Faster, but less
        #  accurate.)
        self.fdSensitivities = False

        # Integrate with log concentrations (to avoid negative concentrations)
        self.integrateWithLogs = False
        
    add_int_times, add_tail_times = True, True
    def full_speed(cls):
        """
        Do not add timepoints to returned trajectories.
        """
        cls.add_int_times, cls.add_tail_times = False, False
    full_speed = classmethod(full_speed)
    def fill_traj(cls):
        """
        Allow the integrator to automatically add timepoints to the interior of
        a trajectory.
        
        This is slower, but may capture the dynamics better.
        """
        cls.add_int_times, cls.add_tail_times = True, False
    fill_traj = classmethod(fill_traj)
    def pretty_plotting(cls):
        """
        Add timepoints to the interior of a trajectory, plus 5% past the last
        timepoint requested.

        This produces the prettiest plots.
        """
        cls.add_int_times, cls.add_tail_times = True, True
    pretty_plotting = classmethod(pretty_plotting)

    #
    # Methods used to build up a network
    #
    def _add_variable(self, var):
        self._checkIdUniqueness(var.id)
        self.variables.set(var.id, var)
        self._makeCrossReferences()

    def add_compartment(self, id, initial_size=1.0, name='', 
                        typical_value=None,
                        is_constant=True, is_optimizable=False):
        """
        Add a compartment to the Network.

        All species must reside within a compartment.
        """
        compartment = Compartment(id, initial_size, name, typical_value, 
                                  is_constant, is_optimizable)
        self._add_variable(compartment)

    def add_species(self, id, compartment, initial_conc=None, 
                    name='', typical_value=None,
                    is_boundary_condition=False, is_constant=False, 
                    is_optimizable=False):
        """
        Add a species to the Network.
        """
        species = Species(id, compartment, initial_conc, name, typical_value,
                          is_boundary_condition, is_constant, is_optimizable)
        self._add_variable(species)

    def add_parameter(self, id, initial_value=1.0, name='',
                      typical_value=None,
                      is_constant=True, is_optimizable=True):
        """
        Add a parameter to the Network.
        """
        parameter = Parameter(id, initial_value, name, is_constant, 
                              typical_value, is_optimizable)
        self._add_variable(parameter)

    def add_event(self, id, trigger, event_assignments={}, delay=0, name=''):
        """
        Add an event to the Network.

        id - id for this event
        trigger - The event firest when trigger passes from False to True.
            Examples: To fire when time becomes greater than 5.0:
                       trigger = 'gt(time, 5.0)'
                      To fire when A becomes less than sin(B/C):
                       trigger = 'lt(A, sin(B/C))'
        event_assignments - A dictionary of assignments to make when the
            event executes.
            Example: To set A to 4.3 and D to B/C
                      event_assignments = {'A': 4.3,
                                           'D': 'B/C'}
        delay - Optionally, assignments may take effect some time after the
            event fires. delay may be a number or math expression
        name - A more detailed name for the event, not restricted to the id
            format
        """
        event = Event(id, trigger, event_assignments, delay, name)
        self._checkIdUniqueness(event.id)
        self.events.set(event.id, event)

    def add_func_def(self, id, variables, math, name=''):
        """
        Add a function definition to the Network.

        id - id for the function definition
        variables - The variables used in the math expression whose
                    values should be subsituted.
        math - The math expression the function definition represents
        name - A more extended name for the definition

        Example:
            To define f(x, y, z) = y**2 - cos(x/z)
            net.add_func_def('my_func', ('x', 'y'), 'y**2 - cos(x/z)')
        """
        func = FunctionDefinitions(id, variables, math, name)
        self._checkIdUniqueness(function.id)
        self.functionDefinitions.set(function.id, function)

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
        self.reactions.set(rxn.id, rxn)

    def add_assignment_rule(self, var_id, rhs):
        """
        Add an assignment rule to the Network.

        A rate rules species that <var_id> = rhs.
        """
        self.assignmentRules.set(var_id, rhs)

    def add_rate_rule(self, var_id, rhs):
        """
        Add a rate rule to the Network.

        A rate rules species that d <var_id>/dt = rhs.
        """
        self.variables.get(var_id).is_constant = False
        self.rateRules.set(var_id, rhs)

    def remove_component(self, id):
        """
        Remove the component with the given id from the Network.

        Components that can be removed are variables, reactions, events,
        function definitions, assignment rules, and rate rules.
        """
        complists = [self.variables, self.reactions, self.functionDefinitions,
                     self.events, self.assignmentRules, self.rateRules]
        for complist in complists:
            # If the id is in a list and has a non-empty name
            if complist.has_key(id):
                complist.remove_by_key(id)

        self._makeCrossReferences()

    def _checkIdUniqueness(self, id):
        """
        Check whether a given id is already in use by this Network.
        """
        if id in self.variables.keys()\
           or id in self.reactions.keys()\
           or id in self.functionDefinitions.keys()\
           or id in self.events.keys()\
           or id == self.id:
            raise ValueError, ('The id %s is already in use!' % id)

    def set_id(self, id):
        """
        Set the id of this Network
        """
        self.id = id

    def get_id(self):
        """
        Get the id of this Network.
        """
        return self.id

    def set_name(self, name):
        """
        Set the name of this Network
        """
        self.id = name 

    def get_name(self):
        """
        Get the name of this Network.
        """
        return self.name

    #
    # Methods to become a 'SloppyCell.Model'
    #

    def Calculate(self, vars, params = None):
        # Add in the times required by all the variables, then convert back to 
        #  a sorted list.
        # Make sure we start from t = 0
        t = sets.Set([0])
        for var, times in vars.items():
            t.union_update(sets.Set(times))
        t = list(t)
        t.sort()

        self.trajectory = self.integrate(t, params, addTimes = True)
    
    def CalculateSensitivity(self, vars, params):
        t = sets.Set([0])

        for var,times in vars.items():
            t.union_update(sets.Set(times))

        t = list(t)
        t.sort()

        self.ddv_dpTrajectory = self.integrateSensitivity(t,params, addTimes = True, rtol = 1.0e-7)
        # we also have the normal trajectory within this trajectory
        indexlength = len(self.dynamicVars) + len(self.assignedVars)
        keyToColumn = KeyedList(zip(self.dynamicVars.keys()
                                    + self.assignedVars.keys(),
                                    range(indexlength)))

        self.trajectory = Trajectory(self, keyToColumn)
        self.trajectory.values = self.ddv_dpTrajectory.values[:,0:indexlength]
        self.trajectory.timepoints = self.ddv_dpTrajectory.timepoints

    def GetName(self):
        return self.id

    def GetParameters(self):
        return KeyedList([(var.id, var.value) for var in
                          self.optimizableVars.values()])

    def GetParameterTypicalValues(self):
        return KeyedList([(var.id, var.typicalValue) for var in
                          self.optimizableVars.values()])

    def GetResult(self, vars):
        result = {}
        times = self.trajectory.timepoints
        for id in vars:
            traj = self.trajectory.getVariableTrajectory(id)
            result[id] = dict(zip(times, traj))

        return result

    def GetSensitivityResult(self, vars):
       opts = self.optimizableVars
       result = {}
       times = self.ddv_dpTrajectory.timepoints
       for id in vars.keys():
           result[id] = {}
           for tIndex,t in enumerate(times) :
               result[id][t] = {}
               for optparams in opts.keys() :
                       result[id][t][optparams] = \
                       self.ddv_dpTrajectory.getVariableTrajectory((id,optparams))[tIndex]
       # note: returns all the timepoints we have, not just
       # the requested ones (which should be a subset of all
       # the timepoints)
       return result

    #
    # The actual integration
    #
    def integrate(self, times, params = None,
                  returnEvents = False, addTimes = True,
                  rtol = None):
        if HAVE_DYNAMICS:
            return Dynamics.integrate(self, times, params)

        self.compile()

        if params is not None:
            self.setOptimizables(params)

        # We can add extra times to try to catch events better
        if addTimes:
            times = sets.Set(times)
            toAdd = scipy.linspace(min(times), max(times), 200)
            times.union_update(sets.Set(toAdd))
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
        if HAVE_DYNAMICS:
            return Dynamics.integrate_sensitivity(self, times, params, rtol)

        self.compile()

        if params is not None:
            self.setOptimizables(params)

        # We can add extra times to try to catch events better
        if addTimes:
            times = sets.Set(times)
            toAdd = scipy.linspace(min(times), max(times), 200)
            times.union_update(sets.Set(toAdd))
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

                stepsize = 1.0e-6*abs(var.typicalValue) + 1.0e-12
                self.set_initial_var_value(id, saved + stepsize)

	        self.resetDynamicVariables()
                tStepped, oaStepped, teStepped, yeStepped, ieStepped = \
                        Integration.Integrate(self, times, rtol = rtol)
                #XXX: Should pull values corresponding to times out of tStepped,
                # oaStepped since param change might mean events fire at 
                # different times or not at all.
		deriv = (oaStepped - oaInitial)/stepsize
		self.set_initial_var_value(id, saved)

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

    def evaluate_expr(self, expr, time=0):
        """
        Evaluate the given expression using the current values of the network
        variables.
        """
        if isinstance(expr, str):
            # We remove beginning and trailing whitespace, just for convenience
            expr = expr.strip()
            expr = self.substituteFunctionDefinitions(expr)
            #expr = self.substituteVariableNames(expr)
            # We substitute float values for all the variable ids in our string
            variables_used = Parsing.extractVariablesFromString(expr)
            variables_used.discard('time')
            while variables_used:
                for id in variables_used:
                    mapping = str(self.variables.getByKey(id).value)
                    expr = Parsing.substituteVariableNamesInString(expr, id, 
                                                                   mapping) 
                variables_used = Parsing.extractVariablesFromString(expr)
                variables_used.discard('time')
            return eval(expr)
        else:
            return expr

    #
    # Methods to get and set object properties
    #
    def set_var_ic(self, id, value):
        """
        Set the initial condition of the variable with the given id.
        """
        if id in self.assignedVars.keys():
            print 'WARNING! Attempt to assign an initial condition to the variable %s, which is determined by an assignment rule. This is a meaningless operation. Instead, change the initial condition of one or more of the components in the rule: %s' % (id, self.assignmentRules.get(id))

        var = self.variables.getByKey(id)
        var.initialValue = value
        if var.is_constant:
            var.value = value
        self.constantVarValues = [var.value for var in
                                  self.constantVars.values()]

    def set_var_value(self, id, value, time=0):
        """
        Set the current stored value of the variable with the given id.
        """
        if id in self.assignedVars.keys():
            print 'WARNING! Attempt to assign a value to the variable %s, which is determined by an assignment rule. This is a meaningless operation. Instead, change the value of one or more of the components in the rule: %s' % (id, self.assignmentRules.get(id))

        var = self.variables.get(id)
        var.value = value
        self.updateAssignedVars(time)

    def setTypicalVariableValue(self, id, value):
        """
        Set the typical value of the variable with the given id.
        """
        self.variables.getByKey(id).typicalValue = value

    def update_optimizable_vars(self, params):
        """
        Update the net's optimizable vars from a passed-in KeyedList, 
         dictionary, or sequence.

        Only those variables that intersect between the net and params are
        changed if a KeyedList or dict is passed in.
        """
        if hasattr(params, 'get'):
            inBoth = sets.Set(self.optimizableVars.keys())
            inBoth = inBoth.intersection(sets.Set(params.keys()))
            for id in inBoth:
                self.set_initial_var_value(id, params.get(id))
        elif len(params) == len(self.optimizableVars):
            for ii, id in enumerate(self.optimizableVars.keys()):
                self.set_initial_var_value(id, params[ii])
        else:
            raise ValueError, 'Passed in parameter set does not have the proper length!'

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

        # We need to update the assigned variables both before and after
        #  we do the pending dynamic vars.
        self.updateAssignedVars(time = 0)

        for id, var in pending:
            self.dynamicVars.getByKey(id).value = \
                    self.evaluate_expr(var.initialValue, 0)
	
	self.updateAssignedVars(time = 0)

    def set_dyn_var_ics(self, values):
        for ii, id in enumerate(self.dynamicVars.keys()):
            self.set_initial_var_value(id, values[ii])

    def updateVariablesFromDynamicVars(self, values, time):
        for ii in range(len(self.dynamicVars)):
            self.dynamicVars[ii].value = values[ii]
        self.updateAssignedVars(time)

    def updateAssignedVars(self, time):
        for id, rhs in self.assignmentRules.items():
            self.assignedVars.getByKey(id).value = self.evaluate_expr(rhs, time)

    def set_var_optimizable(self, id, is_optimizable):
        self.variables.get(id).is_optimizable = is_optimizable

    def set_var_constant(self, id, is_constant):
        self.variables.get(id).is_constant = is_constant

    #
    # Generate the differential equations and functions to calculate them.
    #

    def _makeDiffEqRHS(self):
        # Start with all right-hand-sides set to '0'
        self.diff_eq_rhs = KeyedList([(id, '0') for id
                                    in self.dynamicVars.keys()])

        for id, rxn in self.reactions.items():
            rateExpr = rxn.kineticLaw

	    for reactantId, dReactant in rxn.stoichiometry.items():
                if self.variables.getByKey(reactantId).is_boundary_condition or\
                   self.variables.getByKey(reactantId).is_constant or\
                   self.assignmentRules.has_key(reactantId):
                    # Variables that are boundary conditions, are constant, or
                    #  are assigned aren't modified by reactions, so we move on
                    #  to the next.
                    continue

                reactantIndex = self.dynamicVars.indexByKey(reactantId)
                if (type(dReactant) == types.IntType
                    or type(dReactant) == types.FloatType)\
                   and dReactant != 0:
                    if dReactant == 1:
                        self.diff_eq_rhs[reactantIndex] += \
                                ' + ( %s )' % rateExpr
                    elif dReactant == -1:
                        self.diff_eq_rhs[reactantIndex] += \
                                ' - ( %s )' % rateExpr
                    else:
                        self.diff_eq_rhs[reactantIndex] += ' + %f * ( %s )'\
                                % (dReactant, rateExpr)

                elif type(dReactant) == types.StringType:
                    self.diff_eq_rhs[reactantIndex] += ' + ( %s ) * ( %s )'\
                            % (dReactant, rateExpr)

        # Strip off the initial '0 +/- '
        for id, rhs in self.diff_eq_rhs.items():
            if len(rhs) > 4:
                if rhs[:4] == '0 + ':
                    self.diff_eq_rhs.set(id, rhs[4:])
                elif rhs[:4] == '0 - ':
                    self.diff_eq_rhs.set(id, rhs[2:])

        # Handle the rate rules
        for id, rule in self.rateRules.items():
            self.diff_eq_rhs.set(id, rule)

    def getRatesForDV(self,dvName,time) :
        connectedReactions = {}
        # this is very messy way to do this but just for the moment...
        for paramname in self.parameters.keys() :
            exec(paramname+'='+self.parameters.getByKey(paramname).value.__repr__())
        for dynamicVarName in self.trajectory.keyToColumn.keys() :
            tr = self.trajectory.getVariableTrajectory(dynamicVarName)
            timepoints = self.trajectory.timepoints
            minval = min(abs(timepoints-time)) # closest value to time in timepoints
            minindex = mlab.find(minval==abs(timepoints-time)) # this is actually an array
            chemval = tr[minindex[0]]
            exec(dynamicVarName+'='+chemval.__repr__())

        for id,rxn in self.reactions.items() :
            if dvName in rxn.stoichiometry.keys() :
                if rxn.stoichiometry[dvName] != 0 : # doesn't matter if no change in dvName
                    rateexpr ='('+rxn.stoichiometry[dvName].__repr__()+')*('+rxn.kineticLaw+')'
                    rateexprval = eval(rateexpr)
                    connectedReactions[id] = [rxn.stoichiometry,rateexpr,rateexprval]
        return connectedReactions

    def _make_get_ddv_dt(self):
        self.ddv_dt = scipy.zeros(len(self.dynamicVars), scipy.Float)

        functionBody = 'def get_ddv_dt(self, dynamicVars, time):\n\t'
        functionBody += 'ddv_dt = self.ddv_dt\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)
        functionBody += '\n\t'

        for ii, (id, var) in enumerate(self.dynamicVars.items()):
            rhs = self.diff_eq_rhs.getByKey(id)
            if rhs != '0':
                rhs = self.substituteFunctionDefinitions(rhs)
                functionBody += '# Total derivative of %s wrt time\n\t' % (id)
                functionBody += 'ddv_dt[%i] = %s\n\t' % (ii, rhs)

        functionBody += '\n\n\treturn ddv_dt'

        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))
        return functionBody

    def _make_get_d2dv_ddvdt(self):
        self.d2dv_ddvdt = scipy.zeros((len(self.dynamicVars),
                                     len(self.dynamicVars)), scipy.Float)

        functionBody = 'def get_d2dv_ddvdt(self, dynamicVars, time):\n\t'
        functionBody += 'd2dv_ddvdt = self.d2dv_ddvdt\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)
        functionBody += '\n\t'

        for rhsIndex, rhsId in enumerate(self.dynamicVars.keys()):
            rhs = self.diff_eq_rhs.getByKey(rhsId)
            rhs = self.substituteFunctionDefinitions(rhs)
            # Take all derivatives of the rhs with respect to other 
            #  dynamic variables
            for wrtIndex, wrtId in enumerate(self.dynamicVars.keys()):
                deriv = self.takeDerivative(rhs, wrtId)
                if deriv != '0':
                    functionBody += '# Partial derivative of Ddv_Dt for %s wrt %s\n\t' % (rhsId, wrtId)
                    functionBody += 'd2dv_ddvdt[%i, %i] = %s\n\t' % \
                            (wrtIndex, rhsIndex, deriv)

        functionBody += '\n\treturn d2dv_ddvdt'

        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))
        return functionBody

    def _make_get_d2dv_dovdt(self):
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
                rhs = self.diff_eq_rhs.getByKey(rhsId)
                rhs = self.substituteFunctionDefinitions(rhs)
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

        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))
        return functionBody


    #
    # Methods involved in events
    #

    def fireEvent(self, event, dynamicVarValues, time):
        self.updateVariablesFromDynamicVars(dynamicVarValues, time)

        delay = self.evaluate_expr(event.delay, time)

        executionTime = round(time + delay, 8)
        event.new_values[executionTime] = {}
        # Calculate the values that will get assigned to the variables
        for id, rhs in event.event_assignments.items():
            event.new_values[executionTime][id] = self.evaluate_expr(rhs, time)

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
                delay = self.evaluate_expr(delay, time)

        executionTime = round(time + delay, 8)
        event.new_values[executionTime] = {}

        ddv_dt = self.get_ddv_dt(values[:nDV], time)
        new_values = copy.copy(values[:nDV])
        for lhsId, rhs in event.event_assignments.items():
            lhsIndex = self.dynamicVars.indexByKey(lhsId)

            # Compute and store the new value of the lhs
            newVal = self.evaluate_expr(rhs, time)
            event.new_values[executionTime][lhsId] = newVal
            new_values[lhsIndex] = newVal

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
                dtrigger_ddvValue[dynId] = self.evaluate_expr(dtrigger_ddv, 
                                                              time)

            dtrigger_dovValue = {}
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                dtrigger_dov = self.takeDerivative(trigger, optId)
                dtrigger_dovValue[optId] = self.evaluate_expr(dtrigger_dov,
                                                              time)

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
            rhs = event.event_assignments.get(lhsId, lhsId)
            # Everything below should work if the rhs is a constant, as long
            #  as we convert it to a string first.
            rhs = str(rhs)
            rhs = self.substituteFunctionDefinitions(rhs)

            drhs_ddvValue = {}
            for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                drhs_ddv = self.takeDerivative(rhs, dynId)
                drhs_ddvValue[dynId] = self.evaluate_expr(drhs_ddv, time)

            drhs_dovValue = {}
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                drhs_dov = self.takeDerivative(rhs, optId)
                drhs_dovValue[optId] = self.evaluate_expr(drhs_dov, time)

            drhs_dtime = self.takeDerivative(rhs, 'time')
            drhs_dtimeValue = self.evaluate_expr(drhs_dtime, time)

            # Calculate the perturbations to the sensitivities
            for ovIndex, optId in enumerate(self.optimizableVars.keys()):
                for dvIndex, dynId in enumerate(self.dynamicVars.keys()):
                    event.new_values[executionTime].setdefault((lhsId, optId), 0)
                    indexInOA = dvIndex + (ovIndex + 1)*nDV
                    #dv_ddv * ddv_dov
                    event.new_values[executionTime][(lhsId, optId)] += \
                            drhs_ddvValue[dynId] * values[indexInOA]
                    #dv_ddv * ddv_dt * dTf_dov
                    event.new_values[executionTime][(lhsId, optId)] += \
                            self.get_ddv_dt(values[:nDV], time)[dvIndex]\
                            * drhs_ddvValue[dynId] * dTf_dov[optId]

                event.new_values[executionTime].setdefault((lhsId, optId), 0)

                #dv_dov
                event.new_values[executionTime][(lhsId, optId)] += drhs_dovValue[optId]

                #dv_dtime
                event.new_values[executionTime][(lhsId, optId)] += \
                        drhs_dtimeValue * dTf_dov[optId]

                # dv_dt * dTf_dov
                # note that dc_dt * dTf_dov cancels out
                event.new_values[executionTime][(lhsId, optId)] += \
                        -self.get_ddv_dt(new_values, time)[lhsIndex] * dTf_dov[optId]

        return delay

    def executeEvent(self, event, dynamicVarValues, time):
        for id, value in event.new_values[round(time, 8)].items():
            if type(id) == types.StringType:
                dynamicVarValues[self.dynamicVars.indexByKey(id)] = value

        return dynamicVarValues

    def executeEventAndUpdateDdv_Dov(self, event, inValues, time):
        values = copy.copy(inValues)
        # Copy out the new values for the dynamic variables
        for id, value in event.new_values[round(time, 8)].items():
            if type(id) == types.TupleType:
                nDV = len(self.dynamicVars)
                dvIndex = self.dynamicVars.indexByKey(id[0])
                ovIndex = self.optimizableVars.indexByKey(id[1])
                values[dvIndex + (ovIndex + 1)*nDV] = value

            else:
                values[self.dynamicVars.indexByKey(id)] = value

        return values


    def _make_get_eventValues(self):
        self._eventValues = scipy.zeros(len(self.complexEvents), scipy.Float)
        self._eventTerminals = scipy.zeros(len(self.complexEvents))
        self._eventDirections = scipy.zeros(len(self.complexEvents))

        functionBody = 'def get_eventValues(self, dynamicVars, time):\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)

        for ii, event in enumerate(self.complexEvents.values()):
            rhs = self.substituteFunctionDefinitions(event.trigger)
            functionBody += 'self._eventValues[%i] = %s\n\t' % (ii, rhs)
            functionBody += 'self._eventTerminals[%i] = %s\n\t' % \
                    (ii, event.is_terminal)
            functionBody += 'self._eventDirections[%i] = 1\n\t' % ii

        functionBody += '\n\treturn self._eventValues, self._eventTerminals, self._eventDirections'

        return functionBody

    def _make_root_func(self):
        self._root_func = scipy.zeros(len(self.events), scipy.Float)

        functionBody = 'def root_func(self, dynamicVars, time):\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)

        for ii, event in enumerate(self.events.values()):
            rhs = self.substituteFunctionDefinitions(event.trigger)
            functionBody += 'self._root_func[%i] = %s\n\t' % (ii, rhs)

        functionBody += '\n\treturn self._root_func'

        return functionBody

    def _make_get_eventDerivs(self):
        self._eventDerivValues = scipy.zeros(len(self.complexEvents),
                                             scipy.Float)
        
        functionBody = 'def get_eventDerivs(self, dynamicVars, ddv_dt, time):\n\t'
        functionBody = self.addAssignmentRulesToFunctionBody(functionBody)

        for ii, event in enumerate(self.complexEvents.values()):
            trigger = self.substituteFunctionDefinitions(event.trigger)
            rhs = ''
            for id in self.diff_eq_rhs.keys():
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

        symbolic.saveDiffs(os.path.join(_TEMP_DIR, 'diff.pickle'))

        return functionBody

    #
    # Internally useful things
    #
    def _makeCrossReferences(self):
        """
        Create the cross-reference lists for the Network.
        """
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
            mapping[var.__class__].set(id, var)
            if var.is_constant:
                self.constantVars.set(id, var)
                if var.is_optimizable:
                    self.optimizableVars.set(id, var)
            elif id in self.assignmentRules.keys():
                self.assignedVars.set(id, var)
            else:
                self.dynamicVars.set(id, var)

        self.constantVarValues = [var.value for var in 
                                  self.constantVars.values()]

        self.complexEvents = KeyedList()
        self.timeTriggeredEvents = KeyedList()
        for id, event in self.events.items():
            if event.timeTriggered:
                self.timeTriggeredEvents.set(id, event)
            else:
                self.complexEvents.set(id, event)

    def _exec_dynamic_func(self, func):
        """
        Create the executable function corresponding to func's functionBody.
        """
        function_body = getattr(self, '%s_functionBody' % func)
        exec function_body
        setattr(self, func, 
                types.MethodType(locals()[func], self, Network))
        if HAVE_PSYCO:
            psyco.bind(getattr(self, func))

    def compile(self):
        """
        Create the dynamically-generated functions for this Network.

        Note that if a model involves many copies of the same network, with 
        differences only in events or initial conditions, it will save time
        to compile the base Network before copying it.
        """
        # If the structure of our network has changed, remake all the dynamic
        #  functions
        # Note that this runs at least once, since _last_structure doesn't
        #  exist beforehand.
        curr_structure = self._get_structure()
        if curr_structure != getattr(self, '_last_structure', None):
            self._makeCrossReferences()
            self._makeDiffEqRHS()
            for func in self._dynamic_structure_funcs:
                exec 'self.%s_functionBody = self._make_%s()' % (func, func)
                self._exec_dynamic_func(func)
            self._last_structure = copy.deepcopy(curr_structure)

        if self.events != getattr(self, '_last_events', None):
            for func in self._dynamic_event_funcs:
                exec 'self.%s_functionBody = self._make_%s()' % (func, func)
                self._exec_dynamic_func(func)
            self._last_events = copy.deepcopy(self.events)

    def dynamic_function_from_file(self, filename):
        """
        Load a dynamic function from a file. 
        
        The filename must be <function_name>.py
        """
        f = file(filename, 'r')
        function_body = f.read()
        f.close()

        basename = os.path.basename(filename)
        func = os.path.splitext(basename)[0]
        setattr(self, '%s_functionBody' % func, function_body)
        self._exec_dynamic_func(func)

    def _get_structure(self):
        """
        Return a tuple representing the structure of the Network.

        The tuple contains the functionDefinitions, reactions, and rules for
        the Network. It also contains information on the constancy and
        optimizability of the variables.
        """
        var_struct = {}
        structure = (self.functionDefinitions, self.reactions, 
                     self.assignmentRules, self.rateRules, var_struct)
        for id, var in self.variables.items():
            var_struct[id] = (var.is_constant, var.is_optimizable)

        return structure

    def output_dynamic_functions(self, directory = _TEMP_DIR):
        """
        Output .py files for this Network's dynamic functions into the given
        directory.
        """
        for func in self._dynamic_structure_funcs + self._dynamic_event_funcs:
            f = file(os.path.join(directory, '%s.py' % func), 'w')
            f.write(getattr(self, '%s_functionBody' % func))
            f.close()

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
        assigned_used = Parsing.extractVariablesFromString(input)
        assigned_used.difference_update(sets.Set([wrt]))
        assigned_used.intersection_update(sets.Set(self.assignedVars.keys()))
        # Do the chain rule for those variables
        for id in assigned_used:
            rule = self.assignmentRules.getByKey(id)
            rule = self.substituteFunctionDefinitions(rule)
            d2 = self.takeDerivative(rule, wrt)
            if d2 != '0':
                d = symbolic.Diff(input, id)
                output += ' + %s *(%s)' % (d, d2)

        return output

    def substituteFunctionDefinitions(self, input):
        # extractFunctionsFromString returns a set of tuples of the form
        #  (function name, # of arguments)
        # _standardFuncDefs is keyed by tuples like this, but
        # functionDefinitions is just keyed by function name
        reversedFunctions = [((fd.id, len(fd.variables)), fd) for fd in\
                             self.functionDefinitions.values()]
        reversedFunctions.reverse()
        functionsUsed = Parsing.extractFunctionsFromString(input)
        for (id, numvars), function in reversedFunctions:
            if (id, numvars) in functionsUsed:
                input = Parsing.substituteFunctionIntoString(input, id, 
                                                             function.variables,
                                                             function.math)
                functionsUsed = Parsing.extractFunctionsFromString(input)

        standardIds = sets.Set(self._standardFuncDefs.keys())
        standardUsed = standardIds.intersection(functionsUsed)
        while len(standardUsed) > 0:
            for tup in standardUsed:
                function = self._standardFuncDefs.getByKey(tup)
                input = Parsing.substituteFunctionIntoString(input, function.id,
                                                             function.variables,
                                                             function.math)

            functionsUsed = Parsing.extractFunctionsFromString(input)
            standardUsed = standardIds.intersection(functionsUsed)


        return input

    def copy(self, new_id=None, new_name=None):
        """
        Return a copy of the given network, with an optional new id.
        """
        new_net = copy.deepcopy(self)
        if new_id is not None:
            new_net.set_id(new_id)
        if new_name is not None:
            new_net.set_name(new_name)

        return new_net

    def __getstate__(self):
        odict = copy.copy(self.__dict__)

        for func in self._dynamic_structure_funcs + self._dynamic_event_funcs:
            odict[func] = None

        return odict

    def __setstate__(self, newdict):
        self.__dict__.update(newdict)
        self._makeCrossReferences()
        # We exec all our functions, so we're ready to go.
        for func in self._dynamic_structure_funcs + self._dynamic_event_funcs:
            try:
                self._exec_dynamic_func(func)
            except AttributeError:
                pass


    def _addStandardFunctionDefinitions(self):
        standardFuncDefs = [
                            ('gt', ['x', 'y'],  'x - y'),
                            ('geq', ['x', 'y'],  'x - y'),
                            ('lt', ['x', 'y'],  'y - x'),
                            ('leq', ['x', 'y'],  'y - x'),
                            ('pow', ['x', 'n'],  'x**n'),
                            ('root', ['n', 'x'],  'x**(1./n)'),
                            ('sqrt', ['x'],  'x**(.5)'),
                            ('cot', ['x'],  '1./tan(x)'),
                            ('arccot', ['x'],  'atan(1./x)'),
                            ('coth', ['x'],  '1./tanh(x)'),
                            ('arccoth', ['x'],  'arctanh(1./x)'),
                            ('csc', ['x'],  '1./sin(x)'),
                            ('arccsc', ['x'],  'asin(1./x)'),
                            ('csch', ['x'],  '1./sinh(x)'),
                            ('arccsch', ['x'],  'arcsinh(1./x)'),
                            ('sec', ['x'],  '1./cos(x)'),
                            ('arcsec', ['x'],  'acos(1./x)'),
                            ('sech', ['x'],  '1./cosh(x)'),
                            ('arcsech', ['x'],  'arccosh(1./x)'),
                            ('log', ['b','x'], 'log(x)/log(b)')
                            ]

        self._standardFuncDefs = KeyedList()
        for id, vars, math in standardFuncDefs:
            function = FunctionDefinition(id, vars, math)
            self._checkIdUniqueness(function.id)
            self._standardFuncDefs.set((function.id, 
                                        len(function.variables)), function)


    def get_component_name(self, id, TeX_form=False):
        """
        Return a components's name if it exists, else just return its id.
        """
        # These are all the things that have names (except the network itself)
        complists = [self.variables, self.reactions, self.functionDefinitions,
                     self.events]
        # If we don't find a name, we'll just use the id
        name = id
        for complist in complists:
            # If the id is in a list and has a non-empty name
            if complist.has_key(id) and complist.get(id).name:
                name = complist.get(id).name
                break

        # We can also check the network's name
        if id == self.id and self.name:
            name = self.name

        if TeX_form:
            # If we've got one underscore in the name, use that to indicate a 
            #  subscript
            if name.count('_') == 1:
                sp = name.split('_')
                name = '%s_{%s}' % (sp[0], sp[1])
            else:
                # TeX can't handle more than one _ in a name, so we substitute
                name = name.replace('_', r'\_')

        return name

    # Deprecated functions below.

    def addCompartment(self, id, size=1.0, name='', 
                       typicalValue=False,isConstant=True, isOptimizable=False):
        self.add_compartment(id = id, initial_size = size, name = name,
                             typical_value = typicalValue, 
                             is_constant = isConstant, 
                             is_optimizable = isOptimizable)
        #raise exceptions.DeprecationWarning('Method addCompartment is deprecated, use add_compartment instead.')

    addVariable = _add_variable
    
    def addSpecies(self, id, compartment, initialConcentration=None,
                   name='', typicalValue=None, is_boundary_condition=False,
                   isConstant=False, isOptimizable=False):
        self.add_species(id = id, compartment = compartment, 
                         initial_conc = initialConcentration, 
                         name = name,
                         typical_value = typicalValue, 
                         is_boundary_condition = is_boundary_condition, 
                         is_constant = isConstant,
                         is_optimizable = isOptimizable)

    def addParameter(self, id, value = 0.0, 
                     typicalValue = None, name = '',
                     isConstant = True, isOptimizable = True):
        self.add_parameter(id = id, initial_value = value, name = name, 
                           typical_value = typicalValue,
                           is_constant = isConstant, 
                           is_optimizable = isOptimizable)

    addRateRule = add_rate_rule

    def dyn_var_fixed_point(self, dv0 = None):
        return Dynamics.dyn_var_fixed_point(self, dv0)
    FindFixedPoint = dyn_var_fixed_point

    set_initial_var_value = set_var_ic
    setInitialVariableValue = set_var_ic

    setOptimizables = update_optimizable_vars

    def TeXOutputToFile(self, filename):
        """Output TeX-formatted network eqns to a file."""
        IO.eqns_TeX_file(self, filename)

    def addEvent(self, id, trigger, eventAssignments, delay=0, name=''):
        self.add_event(id = id, trigger = trigger, 
                       event_assignments = eventAssignments, name = name,
                       delay = delay)

    addFunctionDefinition = add_func_def
    addAssignmentRule = add_assignment_rule
