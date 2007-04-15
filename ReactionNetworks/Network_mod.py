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

import logging
logger = logging.getLogger('ReactionNetworks.Network_mod')

import scipy

import SloppyCell
import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList

import SloppyCell.ExprManip as ExprManip
# We load a dictionary of previously-taken derivatives for efficiency
ExprManip.load_derivs(os.path.join(SloppyCell._TEMP_DIR, 'diff.pickle'))
# This will save the diffs dictionary upon exit from the python interpreter
if SloppyCell.my_rank == 0:
    import atexit
    atexit.register(ExprManip.save_derivs, os.path.join(SloppyCell._TEMP_DIR, 
                                                        'diff.pickle'))

import Reactions
import SloppyCell.Collections as Collections

from Components import *
import Dynamics
import Trajectory_mod

import PerfectData

# Optional since it's x86-only.
try:
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

class Network:
    # Here are all the functions we dynamically generate. To add a new one,
    #  called, for example, fooFunc, you need to define _make_fooFunc, which
    #  generates the python code for the function in fooFunc_functionBody.
    _dynamic_structure_funcs = ['res_function',
                                'alg_deriv_func',
                                'dres_dcdot_function', 'dres_dc_function',
                                'dres_dp_function']
    _dynamic_event_funcs = ['root_func']
    _dynamic_funcs = _dynamic_structure_funcs + _dynamic_event_funcs

    
    # These are-predefined functions we want in our working namespace
    def two_arg_log(x, base=scipy.e):
        return scipy.log(x)/scipy.log(base)
    _common_namespace = {'log': two_arg_log,
                         'log10': scipy.log10,
                         'exp': scipy.exp,
                         'cos': scipy.cos,
                         'sin': scipy.sin,
                         'tan': scipy.tan,
                         'acos': scipy.arccos,
                         'asin': scipy.arcsin,
                         'atan': scipy.arctan,
                         'cosh': scipy.cosh,
                         'sinh': scipy.sinh,
                         'tanh': scipy.tanh,
                         'arccosh': scipy.arccosh,
                         'arcsinh': scipy.arcsinh,
                         'arctanh': scipy.arctanh,
                         'exponentiale': scipy.e,
                         'pi': scipy.pi,
                         'scipy': scipy
                         }
    # These are functions we need to create but that should be used commonly
    _standard_func_defs = [('pow', ['x', 'n'],  'x**n'),
                           ('root', ['n', 'x'],  'x**(1/n)'),
                           ('sqrt', ['x'],  'x**(.5)'),
                           ('cot', ['x'],  '1/tan(x)'),
                           ('arccot', ['x'],  'atan(1/x)'),
                           ('coth', ['x'],  '1/tanh(x)'),
                           ('arccoth', ['x'],  'arctanh(1/x)'),
                           ('csc', ['x'],  '1/sin(x)'),
                           ('arccsc', ['x'],  'asin(1/x)'),
                           ('csch', ['x'],  '1/sinh(x)'),
                           ('arccsch', ['x'],  'arcsinh(1/x)'),
                           ('sec', ['x'],  '1/cos(x)'),
                           ('arcsec', ['x'],  'acos(1/x)'),
                           ('sech', ['x'],  '1/cosh(x)'),
                           ('arcsech', ['x'],  'arccosh(1/x)'),
                           ]
    _logical_comp_func_defs = [('gt', ['x', 'y'],  'x > y'),
                               ('geq', ['x', 'y'],  'x >= y'),
                               ('lt', ['x', 'y'],  'x < y'),
                               ('leq', ['x', 'y'],  'x <= y'),
                               ('eq', ['x', 'y'],  'y == x'),
                               ('and_func', ['x', 'y'], 'x and y'),
                               ('or_func', ['x', 'y'], 'x or y'),
                               ('not_func', ['x'], 'not x'),
                               ]
                           
    # Add our common function definitions to a func_strs list.
    _common_func_strs = []
    for id, vars, math in _standard_func_defs:
        var_str = ','.join(vars)
        func = 'lambda %s: %s' % (var_str, math)
        _common_func_strs.append((id, func))
        for ii, wrt in enumerate(vars):
            deriv_id = '%s_%i' % (id, ii)
            func = 'lambda %s: %s' % (var_str, ExprManip.diff_expr(math, wrt))
            _common_func_strs.append((deriv_id, func))
    # Also do the logical functions
    for id, vars, math in _logical_comp_func_defs:
        var_str = ','.join(vars)
        func = 'lambda %s: %s' % (var_str, math)
        _common_func_strs.append((id, func))

    def __init__(self, id, name=''):
        self.id, self.name = id, name

        self.functionDefinitions = KeyedList()
        self.reactions = KeyedList()
        self.assignmentRules = KeyedList()
        self.rateRules = KeyedList()
        self.algebraicRules = KeyedList()
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
        self.algebraicVars = KeyedList()
        self.compartments = KeyedList()
        self.parameters = KeyedList()
        self.species = KeyedList()

        # These are also cross references. complexEvents for events whose
        #  firing depends on variable values, and timeTriggeredEvents for
        #  events that fire at a give time.
        self.complexEvents = KeyedList()
        self.timeTriggeredEvents = KeyedList()

        # All expressions are evaluated, and functions exec'd in self.namespace.
        #  (A dictionary mapping names to objects they represent)
        self.namespace = copy.copy(self._common_namespace)

        # These are the strings we eval to create the functions our Network
        #  defines. We'll evaluate them all to start with, to get the common
        #  ones in there.
        self._func_strs = copy.copy(self._common_func_strs)
        for func_id, func_str in self._func_strs:
            self.namespace[func_id] = eval(func_str, self.namespace, {})

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

    def add_species(self, id, compartment, initial_conc=0, 
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
        event_assignments - A dictionary or KeyedList of assignments to make
                      when the event executes.
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
            net.add_func_def('my_func', ('x', 'y', 'z'), 'y**2 - cos(x/z)')
        """
        func = FunctionDefinition(id, variables, math, name)
        self._checkIdUniqueness(func.id)
        self.functionDefinitions.set(func.id, func)

        # Add the function and its partial derivatives to func_strs
        # Also do the evaluation
        var_str = ','.join(variables)
        func_str = 'lambda %s: %s' % (var_str, math)
        self._func_strs.append((id, func_str))
        self.namespace[id] = eval(func_str, self.namespace)
        for ii, wrt in enumerate(variables):
            diff_id = '%s_%i' % (id, ii)
            func_str = 'lambda %s: %s' % (var_str, ExprManip.diff_expr(math, 
                                                                       wrt))
            self._func_strs.append((diff_id, func_str))
            self.namespace[diff_id] = eval(func_str, self.namespace)

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
        self.set_var_constant(var_id, False)
        self.assignmentRules.set(var_id, rhs)
        self._makeCrossReferences()

    def add_rate_rule(self, var_id, rhs):
        """
        Add a rate rule to the Network.

        A rate rules species that d <var_id>/dt = rhs.
        """
        self.set_var_constant(var_id, False)
        self.rateRules.set(var_id, rhs)
        self._makeCrossReferences()

    def add_algebraic_rule(self, rhs):
        """
        Add an algebraic rule to the Network.

        An algebraic rule specifies that 0 = rhs.
        """
        self.algebraicRules.set(rhs, rhs)
        self._makeCrossReferences()

    def remove_component(self, id):
        """
        Remove the component with the given id from the Network.

        Components that can be removed are variables, reactions, events,
        function definitions, assignment rules, and rate rules.
        """
        complists = [self.variables, self.reactions, self.functionDefinitions,
                     self.events, self.assignmentRules, self.rateRules,
                     self.algebraicRules]
        for complist in complists:
            # If the id is in a list and has a non-empty name
            if complist.has_key(id):
                complist.remove_by_key(id)

        self._makeCrossReferences()

    def _checkIdUniqueness(self, id):
        """
        Check whether a given id is already in use by this Network.
        """
        if id == 'time':
            logger.warn("Specifying 'time' as a variable is dangerous! Are you "
                        "sure you know what you're doing?")
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
    def calculate(self, vars, params=None):
        self.Calculate(vars, params)
        return self.GetResult(vars)

    def Calculate(self, vars, params = None):
        # Add in the times required by all the variables, then convert back to 
        #  a sorted list.
        # Make sure we start from t = 0
        t = sets.Set([0])
        ret_full_traj=False
        for var, times in vars.items():
            if var == 'full trajectory':
                ret_full_traj = True
                t.union_update(sets.Set(times))
            elif self.variables.has_key(var):
                t.union_update(sets.Set(times))
            else:
                raise ValueError('Unknown variable %s requested from network %s'
                                 % (var, self.get_id()))

        # This takes care of the normal data points
        t = list(t)
        t.sort()
        self.trajectory = self.integrate(t, params=params, addTimes = ret_full_traj)
        
    def CalculateSensitivity(self, vars, params):
        t = sets.Set([0])

        for var,times in vars.items():
            t.union_update(sets.Set(times))

        t = list(t)
        t.sort()

        self.ddv_dpTrajectory = self.integrateSensitivity(t,params=params, 
                                                          addTimes=True)
        self.trajectory = self.ddv_dpTrajectory

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
            if id == 'full trajectory':
                self.trajectory.build_interpolated_traj()
                result[id] = self.trajectory
            elif self.variables.has_key(id):
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

    def integrate(self, times, params = None,
                  returnEvents = False, addTimes = True,
                  rtol = None):
        if self.add_tail_times or addTimes:
            times = scipy.concatenate((times, [1.05*times[-1]]))
        return Dynamics.integrate(self, times, params=params, 
                                  fill_traj=(self.add_int_times or addTimes))

    def integrateSensitivity(self, times, params = None,
                             returnEvents = False, addTimes = True,
                             rtol=1e-6):
        if self.add_tail_times:
            times = scipy.concatenate((times, [1.05*times[-1]]))
        return Dynamics.integrate_sensitivity(self, times, params, rtol,
                                              fill_traj=self.add_int_times)

    def _sub_for_piecewise(self, expr, time):
        """
        Runs through expr and replaces all piecewise expressions by the
        clause that is currently active.
        """
        funcs_used = ExprManip.extract_funcs(expr)
        for func, args in funcs_used:
            if func == 'piecewise':
                ast = ExprManip.strip_parse(expr)
                ast = self._sub_for_piecewise_ast(ast, time)
                return ExprManip.ast2str(ast)
        else:
            return expr
    
    def _sub_for_piecewise_ast(self, ast, time):
        if isinstance(ast, ExprManip.AST.CallFunc)\
           and ExprManip.ast2str(ast.node) == 'piecewise':
            # If our ast is a piecewise function
            conditions = [cond for cond in ast.args[1::2]]
            clauses = [clause for clause in ast.args[:-1:2]]
            if len(ast.args)%2 == 1:
                # odd # of arguments implies otherwise clause
                otherwise = ast.args[-1]
            else:
                otherwise = None
        
            for cond_ast, expr_ast in zip(conditions, clauses):
                # For each condition, check whether or not it is True.
                # If so, return the corresponding clause.
                subbed_cond = self._sub_for_piecewise_ast(cond_ast, time=time)
                cond_str = ExprManip.ast2str(subbed_cond)
                if self.evaluate_expr(cond_str, time):
                    return self._sub_for_piecewise_ast(expr_ast, time=time)
            else:
                # If none of our conditions were True, return the otherwise
                # clause if we have one.
                if otherwise is not None:
                    return self._sub_for_piecewise_ast(otherwise, time=time)
                else:
                    raise ValueError('No True condition and no otherwise '
                                     " clause in '%s'."
                                     % ExprManip.ast2str(ast))
        return ExprManip.AST.recurse_down_tree(ast, self._sub_for_piecewise_ast,
                                           (time,))

    def evaluate_expr(self, expr, time=0):
        """
        Evaluate the given expression using the current values of the network
        variables.
        """
        if isinstance(expr, str):
            # We create a local_namespace to evaluate the expression in that
            #  maps variable ids to their current values
            expr = self._sub_for_piecewise(expr, time)
            vars_used = ExprManip.extract_vars(expr)
            var_vals = [(id, self.get_var_val(id)) for id in vars_used
                        if id != 'time']
            local_namespace = dict(var_vals)
            local_namespace['time'] = time
            local_namespace.update(self.namespace)
            # We strip whitespace, just for convenience
            return eval(expr.strip(), local_namespace, {})
        else:
            return expr

    #
    # Methods to get and set object properties
    #
    def _get_var_attr(self, id, attr):
        if self.variables.has_key(id):
            return getattr(self.get_variable(id), attr)
        else:
            raise ValueError, 'Id %s not found in network.' % id

    def set_var_typical_val(self, id, value):
        """
        Set the typical value for a variable.
        """
        var = self.get_variable(id)
        var.typicalValue = value

    def get_var_typical_val(self, id):
        """
        Return the typical value for a variable.
        """
        return self._get_var_attr(id, 'typicalValue')

    def get_var_typical_vals(self):
        """
        Return the variable typical values as a KeyedList.
        """
        return KeyedList([(id, self.get_var_typical_val(id)) for id in
                          self.variables.keys()])

    def set_var_ic(self, id, value, warn=True):
        """
        Set the initial condition of the variable with the given id.
        """
        if warn and id in self.assignedVars.keys():
            logger.warn('WARNING! Attempt to assign an initial condition to '
                        'the variable %s, which is determined by an assignment '
                        'rule. This is a meaningless operation. Instead, '
                        'change the initial condition of one or more of the '
                        'components in the rule: %s' % 
                        (id, self.assignmentRules.get(id)))

        var = self.get_variable(id)
        var.initialValue = value
        if var.is_constant:
            var.value = value
        self.constantVarValues = [self.evaluate_expr(var.value) for var in
                                  self.constantVars.values()]
        self.constantVarValues = scipy.array(self.constantVarValues)

    def get_var_ic(self, id):
        """
        Return the initial condition for a variable
        """
        return self._get_var_attr(id, 'initialValue')

    def get_var_ics(self):
        """
        Return the variable initial conditions as a KeyedList
        """
        return KeyedList([(id, self.get_var_ic(id)) for id in
                          self.variables.keys()])

    def set_var_vals(self, kl, time = 0):
        """
        Set current variable values from a KeyedList or dictionary.
        """
        for id, value in kl.items():
            self.set_var_val(id, value, time, warn=False, do_assignments=False)
        self.updateAssignedVars(time)

    def get_var_val(self, id):
        """
        Return the current value of a variable
        """
        val = self._get_var_attr(id, 'value')
        return self.evaluate_expr(val)

    def get_var_vals(self):
        """
        Return the current variable values as a KeyedList
        """
        return KeyedList([(id, self.get_var_val(id)) for id in
                          self.variables.keys()])

    def get_initial_velocities(self) :
        """
        Returns the vector field evaluated at the initial conditions
        """
        ics = [self.evaluate_expr(self.getInitialVariableValue(dvid))
               for dvid in self.dynamicVars.keys()]
        # the evaluate_expr is in case an initial condition is a parameter
        initialv = Dynamics.find_ics(ics, 0*ics + 1, 0, self.res_function,
                                     self.alg_deriv_func,
                                     self._dynamic_var_algebraic, 1e-6)
        return initialv

    def set_var_ics(self, kl):
        """
        Set variable initial conditions from a KeyedList or dictionary.
        """
        for id, value in kl.items():
            if id in self.variables.keys():
                self.set_var_ic(id, value, warn=False)

    def set_var_val(self, id, val, time=0, warn=True, do_assignments=True):
        """
        Set the current stored value of the variable with the given id.
        """

        if warn and self.assignedVars.has_key(id):
            logger.warn('WARNING! Attempt to assign a value to the variable '
                        '%s, which is determined by an assignment rule. This '
                        'is a meaningless operation. Instead, change the value '
                        'of one or more of the components in the rule: %s'
                        % (id, self.assignmentRules.get(id)))

        var = self.get_variable(id)
        var.value = val
        if do_assignments:
            self.updateAssignedVars(time)

    def set_var_typical_val(self, id, val):
        """
        Set the typical value of the variable with the given id.
        """
        self.get_variable(id).typicalValue = val

    setTypicalVariableValue = set_var_typical_val

    def get_variable(self, id):
        """
        Return the class instance with a given variable id.
        """
        var = self.variables.get(id)
        if var:
            return var
        else:
            raise KeyError('Variable %s not found in network %s!'
                           % (id, self.get_id()))

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
            raise ValueError('Passed in parameter set does not have the proper '
                             'length!')

    getInitialVariableValue = get_var_ic

    def getDynamicVarValues(self):
        # We need to evaluate_expr here to handle ones that are assigned to
        #  parameters
        return scipy.array([self.evaluate_expr(var.value) 
                            for var in self.dynamicVars.values()])

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
            if hasattr(values, 'get'):
                self.set_var_ic(id, values.get(id))
            else:
                self.set_var_ic(id, values[ii])

    def updateVariablesFromDynamicVars(self, values, time):
        for ii in range(len(self.dynamicVars)):
            self.dynamicVars[ii].value = values[ii]
        self.updateAssignedVars(time)

    def updateAssignedVars(self, time):
        for id, rhs in self.assignmentRules.items():
            self.assignedVars.getByKey(id).value = self.evaluate_expr(rhs, time)

    def set_var_optimizable(self, id, is_optimizable):
        self.get_variable(id).is_optimizable = is_optimizable
        self._makeCrossReferences()

    def set_var_constant(self, id, is_constant):
        self.get_variable(id).is_constant = is_constant
        self._makeCrossReferences()

    def get_var_constant(self, id):
        return self.get_variable(id).is_constant

    #
    # Generate the differential equations and functions to calculate them.
    #

    def _makeDiffEqRHS(self):
        logger.debug('Making diff equation rhs')
        diff_eq_terms = {}

        for rxn_id, rxn in self.reactions.items():
            rateExpr = rxn.kineticLaw
            logger.debug('Parsing reaction %s.' % rxn_id)

	    for reactantId, dReactant in rxn.stoichiometry.items():
                if self.get_variable(reactantId).is_boundary_condition or\
                   self.get_variable(reactantId).is_constant or\
                   self.assignmentRules.has_key(reactantId):
                    # Variables that are boundary conditions, are constant, or
                    #  are assigned aren't modified by reactions, so we move on
                    #  to the next.
                    continue

                term = '(%s) * (%s)'  % (dReactant, rateExpr)
                diff_eq_terms.setdefault(reactantId, [])
                diff_eq_terms[reactantId].append(term)

        self.diff_eq_rhs = KeyedList()
        for id in self.dynamicVars.keys():
            if self.rateRules.has_key(id):
                self.diff_eq_rhs.set(id, self.rateRules.get(id))
            else:
                # We use .get to return a default of ['0']
                rhs = '+'.join(diff_eq_terms.get(id, ['0']))
                self.diff_eq_rhs.set(id, ExprManip.simplify_expr(rhs))

    def _make_alg_deriv_func(self):
        # This function is used when we want to calculate consistent derivatives
        #  for algebraic variables. It returns the derivative wrt time of
        #  all the algebraic rules. Since the rules are <0=rule>, these should
        #  all be zero when we have consistent values for the variable derivs.
        body = []
        body.append('def alg_deriv_func(self, alg_yp, dynamicVars, yp, time):')
        body.append('alg_derivs = scipy.empty(%i, scipy.float_)'
                    % len(self.algebraicRules))
        self._add_assignments_to_function_body(body)
        body.append('')

        for rule_ii, rule in enumerate(self.algebraicRules.values()):
            # Each rhs is a sum of drule/dA * dA/dt + drule/dB * dB/dt + ...
            rhs_terms = []
            rhs_terms.append(self.takeDerivative(rule, 'time'))
            for dyn_id in self.dynamicVars.keys():
                deriv = self.takeDerivative(rule, dyn_id)
                if deriv != '0':
                    if self.algebraicVars.has_key(dyn_id):
                        index = self.algebraicVars.index_by_key(dyn_id)
                        rhs_terms.append('(%s)*alg_yp.item(%i)' % (deriv,index))
                    else:
                        index = self.dynamicVars.index_by_key(dyn_id)
                        rhs_terms.append('(%s)*yp.item(%i)' % (deriv, index))
            rhs = ' + '.join(rhs_terms)
            body.append('alg_derivs[%i] = %s' % (rule_ii, rhs))

        body.append('')
        body.append('return alg_derivs')
        return '\n    '.join(body)


    def _make_res_function(self):
        self.residual = scipy.zeros(len(self.dynamicVars), scipy.float_)
        body = []
        body.append('def res_function(self, time, dynamicVars, yprime, ires):')
        body.append('residual = self.residual')
        self._add_assignments_to_function_body(body)
        body.append('')

        # Make a list of algebraic rules for accessing later
        algebraicRuleList = self.algebraicRules.values()

        # We keep a list of terms in our residual function for use building the 
        #  various derivative functions.
        self._residual_terms = []
        # This list tells us whether a give dynamic variable is algebraic or not
        self._dynamic_var_algebraic = []
        # Loop over everything in the dynamicVars list
        for ii, (id, var) in enumerate(self.dynamicVars.items()):
            if self.algebraicVars.has_key(id):
                # It's an algebraic equation. Pop an algebraic rule off the
                # list.
                rhs = algebraicRuleList.pop()
                body.append('# Residual function corresponding to an '
                            'algebraic  equation')
                body.append('residual[%i] = %s' % (ii, rhs))
                self._residual_terms.append((rhs, None))
                self._dynamic_var_algebraic.append(-1)
            else:
                rhs = self.diff_eq_rhs.getByKey(id)
                self._dynamic_var_algebraic.append(+1)
                body.append('# Residual function for %s' % id)
                if rhs != '0':
                    body.append('residual[%i] = %s - yprime[%i]' % 
                                (ii, rhs, ii))
                    self._residual_terms.append((rhs, id))
                else:
                    body.append('residual[%i] = - yprime[%i]' % (ii, ii))
                    self._residual_terms.append((None, id))

        body.append('')
        body.append('return residual')
        return '\n\t'.join(body)

    def _make_dres_dcdot_function(self):
        self.dres_dcdot = scipy.zeros((len(self.residual),
                                       len(self.residual)), 
                                      scipy.float_)
        body = []
        body.append('def dres_dcdot_function(self, time, dynamicVars, yprime, '
                    'indices = None):')
        body.append('dres_dcdot = self.dres_dcdot')
        self._add_assignments_to_function_body(body)
        body.append('')

        for wrt_ii, wrt in enumerate(self.dynamicVars.keys()):
            body.append('if (indices is None) or (%i in indices):' % wrt_ii)
            for res_ii, (rhs, dt) in enumerate(self._residual_terms):
                if dt == wrt:
                    body.append('\t# Derivative of residual term %i wrt %s_dot'
                                % (res_ii, wrt))
                    body.append('\tdres_dcdot[%i, %i] = -1' % 
                                (res_ii, wrt_ii))
            body.append('\tpass')

        body.append('return dres_dcdot')

        return '\n\t'.join(body)

    def _make_dres_dc_function(self):
        self.dres_dc = scipy.zeros((len(self.residual), len(self.residual)), 
                                   scipy.float_)
        body = []
        body.append('def dres_dc_function(self, time, dynamicVars, yprime, '
                    'indices = None):')
        body.append('dres_dc = self.dres_dc')
        self._add_assignments_to_function_body(body)
        body.append('')

        for wrt_ii, wrt in enumerate(self.dynamicVars.keys()):
            body.append('if (indices is None) or (%i in indices):' % wrt_ii)
            for res_ii, (rhs, dt) in enumerate(self._residual_terms):
                if rhs is None:
                    continue
                deriv = self.takeDerivative(rhs, wrt)
                if deriv != '0':
                    if dt is None:
                        # It was an algebraic rule.
                        body.append('\t# Algebraic equation %s wrt %s' %
                                    (rhs, wrt))
                    else:
                        body.append('\t# Residual function for %s wrt %s' %
                                    (dt, wrt))
                    body.append('\tdres_dc[%i, %i] = %s' % 
                                (res_ii, wrt_ii, deriv))
            body.append('\tpass')

        body.append('return dres_dc')

        return '\n\t'.join(body)

    def _make_dres_dp_function(self):
        self.dres_dp = scipy.zeros((len(self.residual),
                                    len(self.optimizableVars)), 
                                   scipy.float_)
        body = []
        body.append('def dres_dp_function(self, time, dynamicVars, yprime, '
                    'indices = None):')
        body.append('dres_dp = self.dres_dp')
        self._add_assignments_to_function_body(body)
        body.append('')

        for wrt_ii, wrt in enumerate(self.optimizableVars.keys()):
            body.append('if (indices is None) or (%i in indices):' % wrt_ii)
            for res_ii, (rhs, dt) in enumerate(self._residual_terms):
                if rhs is None:
                    continue
                deriv = self.takeDerivative(rhs, wrt)
                if deriv != '0':
                    if dt is None:
                        # It was an algebraic rule.
                        body.append('\t# Algebraic equation %s wrt %s' %
                                    (rhs, wrt))
                    else:
                        body.append('\t# Residual function for %s wrt %s' %
                                    (dt, wrt))
                    body.append('\tdres_dp[%i, %i] = %s' % 
                                (res_ii, wrt_ii, deriv))
            body.append('\tpass')

        body.append('return dres_dp')
        return '\n\t'.join(body)

    def _add_assignments_to_function_body(self, body):
        """
        Adds the assignment rules for this Network to the list of lines that
        are in body.
        """
        body.append('constants = self.constantVarValues')

        # We loop to assign our constantVarValues and our dynamicVars
        #  to local variables for speed (avoid repeated accesses) and
        #  for readability
        for arg, var_names in zip(['constants', 'dynamicVars'],
                                  [self.constantVars.keys(),
                                   self.dynamicVars.keys()]):
            # We protect ourselves with a 'try, except' clause
            #  in case passed in values are not an array, but e.g. a list
            #  or KeyedList.
            body.append('try:')
            for ii, id in enumerate(var_names):
                body.append('\t%s = %s.item(%i)' % (id, arg, ii))
            body.append('except (AttributeError, TypeError):')
            for ii, id in enumerate(var_names):
                body.append('\t%s = %s[%i]' % (id, arg, ii))
            body.append('')

        for variable, math in self.assignmentRules.items():
            body.append('%s = %s' % (variable, math))

        return body

    #
    # Methods involved in events
    #

    def fireEvent(self, event, dynamicVarValues, time):
        self.updateVariablesFromDynamicVars(dynamicVarValues, time)
        delay = self.evaluate_expr(event.delay, time)

        return delay

    def executeEvent(self, event, time_fired, y_fired, y_current, time_current):
        # The values assigned depend on the state of the network at the time
        #  the event was >fired<. We go back to that time here...
        self.updateVariablesFromDynamicVars(y_fired, time_fired)
        new_values = {}
        # Calculate the values that will get assigned to the variables
        for id, rhs in event.event_assignments.items():
            new_values[id] = self.evaluate_expr(rhs, time_fired)

        # Make all the relevant assignments
        for id, value in new_values.items():
            y_current[self.dynamicVars.indexByKey(id)] = value
        # Update our network with the new values, just for consistency.
        self.updateVariablesFromDynamicVars(y_current, time_current)

        return y_current

    def _smooth_trigger(self, trigger):
        ast = ExprManip.strip_parse(trigger)
        if len(ast.ops) != 1:
            raise ValueError('Comparison has more than one operation in '
                             'clause %s!' % trigger)
        lhs = ExprManip.ast2str(ast.expr)
        rhs = ExprManip.ast2str(ast.ops[0][1])
        if ast.ops[0][0] in ['>', '>=']:
            return '%s - %s' % (lhs, rhs)
        elif ast.ops[0][0] in ['<', '<=']:
            return '-%s + %s' % (lhs, rhs)
        else:
            raise ValueError('Comparison %s in triggering clause is not '
                             '<, <=, >, or >=!')

    def executeEventAndUpdateSens(self, holder, ysens_pre_exec, opt_var):
        # We extract all our relevant quantities from the event_info holder
        event = holder.event
        time_fired = holder.time_fired
        ysens_fired = holder.ysens_fired
        yp_fired = holder.yp_fired
        clause_index = holder.clause_index
        time_exec = holder.time_exec
        yp_pre_exec = holder.yp_pre_exec
        y_post_exec = holder.y_post_exec
        yp_post_exec = holder.yp_post_exec

        # This is subtle... If this event is the result of a chain, then the
        #  firing time for the event is determined not by its trigger, but
        #  by the execution time of the event that caused it to fire.
        # We walk up the list of chained events, so that finally we set the
        #  d(time of firing) for this event to d(time of execution) for the
        #  start of that chain of events.
        dtf_dp = None
        link = holder
        while link.chained_off_of is not None:
            dtf_dp = link.chained_off_of.dte_dp
            link = link.chained_off_of
        
        ysens_post_exec = copy.copy(ysens_pre_exec)
        # Fill in the new values for the normal y variables
        N_dv = len(self.dynamicVars)
        ysens_post_exec[:N_dv] = y_post_exec

        # We set our network to the state at which the event fired, so we
        #  can use evaluate_expr for the calculations of d(firing time)
        self.updateVariablesFromDynamicVars(ysens_fired, time_fired)
        # We need to do this, if we don't have a chained event.
        if dtf_dp is None:
            # We convert our trigger to smooth function we can take the 
            # derivative of so we can calculate sensitivity of firing time.
            # The trigger is passed to us because it might be a clause of a
            # more complicated trigger.
            trigger = self.event_clauses[clause_index]
            trigger = self._smooth_trigger(trigger)

            # Compute the sensitivity of our firing time
            #
            # numerator = sum_dv dtrigger_ddv ddv_dp + dtrigger_dp
            # denominator = sum_dv dtrigger_ddv ddv_dt + dtrigger_dt
            # dtf_dp = -(numerator)/(denominator)
            numerator = 0
            denominator = 0
            # First term in numerator: sum_dv dtrigger_ddv ddv_dp 
            # First term in denominator: sum_dv dtrigger_ddv ddv_dt 
            for ii, id in enumerate(self.dynamicVars.keys()):
                dtrigger_ddv = self.takeDerivative(trigger, id)
                dtrigger_ddv = self.evaluate_expr(dtrigger_ddv, time_fired)
                # This is the index of ddv_dp in the y array
                index_in_y = ii + N_dv
                ddv_dp = ysens_fired[index_in_y]
                ddv_dt = yp_fired[ii]
                numerator += dtrigger_ddv * ddv_dp
                denominator += dtrigger_ddv * ddv_dt

            # Second term in numerator: dtrigger_dp
            dtrigger_dp = self.takeDerivative(trigger, opt_var)
            dtrigger_dp = self.evaluate_expr(dtrigger_dp, time_fired)
            numerator += dtrigger_dp

            # Second term in denominator: dtrigger_dt
            dtrigger_dt = self.takeDerivative(trigger, 'time')
            dtrigger_dt = self.evaluate_expr(dtrigger_dt, time_fired)
            denominator += dtrigger_dt

            dtf_dp = -numerator/denominator

        # Now compute the sensitivity of the delay
        #
        # ddelay/dp = ddelay/dp + ddelay_dy*dy/dp + ddelay_dy*dy/dt*dtf/dp
        #             + ddelay/dt*dtf/dp
        delay = str(event.delay)

        # First term
        ddelay_dp = self.takeDerivative(delay, opt_var)
        ddelay_dp = self.evaluate_expr(ddelay_dp, time_fired)

        # Second and third terms: ddelay_dy*dy/dp + ddelay_dy*dy/dt*dtf/dp
        for ii, id in enumerate(self.dynamicVars.keys()):
            ddelay_dy = self.takeDerivative(delay, id)
            ddelay_dy = self.evaluate_expr(ddelay_dy, time_fired)
            # This is the index of ddv_dp in the y array
            index_in_y = ii + N_dv
            dy_dp = ysens_fired[index_in_y]
            dy_dt = yp_fired[ii]
            ddelay_dp += ddelay_dy * dy_dp
            ddelay_dp += ddelay_dy * dy_dt * dtf_dp

        # Fourth term: ddelay/dt*dtf/dp
        ddelay_dt = self.takeDerivative(delay, 'time')
        ddelay_dt = self.evaluate_expr(ddelay_dt, time_fired)
        ddelay_dp += ddelay_dt * dtf_dp

        # This is the delay in the time of execution
        dte_dp = dtf_dp + ddelay_dp

        # We store d(time of execution) for use by other events that may be
        #  chained to this one.
        holder.dte_dp = dte_dp
        
        # Now compute the sensitivity of each of our new values
        for y_ii, y_id in enumerate(self.dynamicVars.keys()):
            # Let's say we shift parameters by dp and correspondingly the firing
            #  time of the events shifts by dtf. Then the total change in
            #  y divided by dp  is
            #  dy/dp is dy/dt(pre_exec)*dte/dp + dbump/dp 
            #           - dy/dt(post_exec)*dte/dp
            # Here bump = assigned value - prior value, and
            #  dy/dt(pre_exec) is dy/dt right before the event >executes<.
            #  dy/dt(post_exec) is dy/dt right after the event >executes<.
            # Also dbump/dp = da/dp + da/dy * dy/dp + da/dy * dy/dt * dtf/dp
            #                   + da/dt * dtf/dp - dy/dp
            
            # First we calculate dbump_dp
            # Note that in the logic I've used above, a is y even if there's
            #  no event assignment as long as I evaluate everything at the time
            #  of execution.
            if y_id in event.event_assignments:
                time_bumped = time_fired
                ysens_bumped = ysens_fired
                yp_bumped = yp_fired
            else:
                time_bumped = time_exec
                ysens_bumped = ysens_pre_exec
                yp_bumped = yp_pre_exec
            # We update our Network so evaluate_expr works.
            self.updateVariablesFromDynamicVars(ysens_bumped, time_bumped)

            # If there's no explicit event assignment, we'll have e.g. 'x' = 'x'
            a = event.event_assignments.get(y_id, y_id)
            a = str(a)
            # We only need to deal with whatever clause in the piecewise is
            #  applicable.
            a = self._sub_for_piecewise(a, time_bumped)

            dbump_dp = 0

            # First term: da/dp 
            da_dp = self.takeDerivative(a, opt_var)
            da_dp = self.evaluate_expr(da_dp, time_bumped)
            dbump_dp += da_dp

            # Second term: da/dy * dy/dp 
            # Third term: da/dy * dy/dt * dtf/dp 
            # = da/dy * (dy/dp + dy/dt * dtf/dp)
            for other_ii, other_id in enumerate(self.dynamicVars.keys()):
                da_dother = self.takeDerivative(a, other_id)
                da_dother = self.evaluate_expr(da_dother, time_bumped)

                # This is the index of dother_dp in the y array
                index_in_y = other_ii + N_dv
                dother_dp = ysens_bumped[index_in_y]

                dother_dt = yp_bumped[other_ii]

                dbump_dp += da_dother * (dother_dp + dother_dt * dtf_dp)

            # Fourth term: da/dt * dtf/dp
            da_dt = self.takeDerivative(a, 'time')
            da_dt = self.evaluate_expr(da_dt, time_bumped)
            dbump_dp += da_dt * dtf_dp

            # Fifth term:  -dy/dp = dy/dt * dtf/dp
            dbump_dp -= yp_bumped[y_ii] * dtf_dp

            dy_dp = dbump_dp
            # Now the dy/dt(pre_exec) * dte/dp term
            dy_dp += yp_pre_exec[y_ii] * dte_dp
            # Finally, -dy/dt(post_exec) * dte/dp
            # We calculated dy/dt(new) up above in this function
            dy_dp -= yp_post_exec[y_ii] * dte_dp

            index_in_y = y_ii + N_dv
            ysens_post_exec[index_in_y] = dy_dp

        return ysens_post_exec
        
    def _make_root_func(self):
        len_root_func = 0

        body = []
        body.append('def root_func(self, dynamicVars, time):')
        self._add_assignments_to_function_body(body)

        self.event_clauses = []
        for ii, event in enumerate(self.events.values()):
            trigger = event.trigger
            for func_name, func_vars, func_expr in self._logical_comp_func_defs:
                trigger = ExprManip.sub_for_func(trigger, func_name, func_vars,
                                                 func_expr)

            body.append('self._root_func[%i] = (%s) - 0.5'
                        % (len_root_func, trigger))
            self.event_clauses.append(trigger)
            len_root_func += 1

        for ii, event in enumerate(self.events.values()):
            event.sub_clause_indices = []
            trigger = event.trigger
            for func_name, func_vars, func_expr in self._logical_comp_func_defs:
                trigger = ExprManip.sub_for_func(trigger, func_name, func_vars,
                                                 func_expr)
            comparisons = ExprManip.extract_comps(trigger)
            if len(comparisons) > 1:
                for comp in comparisons:
                    body.append('self._root_func[%i] = (%s) - 0.5\n\t' 
                                % (len_root_func, comp))
                    self.event_clauses.append(comp)
                    event.sub_clause_indices.append(len_root_func)
                    len_root_func += 1

        self._root_func = scipy.zeros(len_root_func, scipy.float_)

        body.append('')
        body.append('return self._root_func')

        return '\n\t'.join(body)

    # ddaskr_root(...) is the function for events for the daskr integrator.
    def ddaskr_root(self, t, y, yprime):
        # note that the order of inputs is reversed in root_func
        return self.root_func(y, t) 

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
        self.algebraicVars = KeyedList()

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

        self.constantVarValues = [self.evaluate_expr(var.value) for var in 
                                  self.constantVars.values()]
        self.constantVarValues = scipy.array(self.constantVarValues)

        self.complexEvents = KeyedList()
        self.timeTriggeredEvents = KeyedList()
        for id, event in self.events.items():
            if event.timeTriggered:
                self.timeTriggeredEvents.set(id, event)
            else:
                self.complexEvents.set(id, event)

        # Collect all variables that are explicitly in algebraic rules
        vars_in_alg_rules = sets.Set()
        for rule in self.algebraicRules:
            vars_in_alg_rules.union_update(ExprManip.extract_vars(rule))

        # Now replace all the assigned variables with the variables they
        # actually depend on. This takes a while loop because assigned
        # vars may be functions of other assigned vars.
        assigned_in_alg =  vars_in_alg_rules.intersection(sets.Set(
                            self.assignedVars.keys()))
        # while there are still assignment variables that we have not
        # expanded in terms of their definitions
        while assigned_in_alg:
            # Replace that assigned_var with the variables it depends on
            # in the running list of algebraic rules
            for assigned_var in assigned_in_alg:
                vars_in_alg_rules.union_update(ExprManip.extract_vars(
                    self.assignmentRules.get(assigned_var)))
                vars_in_alg_rules.remove(assigned_var)
            # update the list of assignment variables that appear in the
            # algebraic rule
            assigned_in_alg =  vars_in_alg_rules.intersection(
                sets.Set(self.assignedVars.keys()))

        # At this point, vars_in_alg_rules should contain all the variables the
        # algebraic rules depend on, including implicit dependencies through
        # assignment rules. Now we filter out all the things we know aren't
        # algebraic vars. First we filter out everything that we already
        # know isn't a dynamic variable.
        vars_in_alg_rules.intersection_update(sets.Set(self.dynamicVars.keys()))

        # remove the reaction variables
        for rxn in self.reactions:
            for chem, value in rxn.stoichiometry.items():
                if value != 0 and (chem in vars_in_alg_rules):
                    vars_in_alg_rules.remove(chem)

        # remove the rate variables
        for var in self.rateRules.keys():
            if (var in vars_in_alg_rules):
                vars_in_alg_rules.remove(var)

        # remove the event variables
        for e in self.events:
            for var in e.event_assignments.keys():
                if (var in vars_in_alg_rules):
                    vars_in_alg_rules.remove(var)

        # Check that we have at most same number of algebraic rules and
        # algebraic vars. Set the algebraicVars list to the list we just
        # compiled
        self.algebraicVars = KeyedList([(id, id) for id
                                        in vars_in_alg_rules])

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
        #  exist beforehand. Also note that __setstate__ will exec all our
        #  dynamic functions upon copying.
        curr_structure = self._get_structure()
        last_structure = getattr(self, '_last_structure', None)
        structure_changed = (curr_structure != last_structure)

        if self.functionDefinitions != getattr(self, '_last_funcDefs', None)\
           or structure_changed:
            logger.debug('Network %s: compiling function defs.' % self.id)
            # Eval our function definitions. Note that these need to be done
            #  in order, and at each point we need to refer back to 
            #  self.namespace, so that, for example, if f(x) = 1  + g(x)
            #  the evaluation of f has the definition of g.
            self.namespace = copy.copy(self._common_namespace)
            for func_id, func_str in self._func_strs:
                self.namespace[func_id] = eval(func_str, self.namespace, {})
            self._last_funcDefs = copy.deepcopy(self.functionDefinitions)

        if curr_structure != getattr(self, '_last_structure', None):
            logger.debug('Network %s: compiling structure.' % self.id)
            self._makeCrossReferences()
            # Check whether system appears well-determined.
            if len(self.algebraicVars) > len(self.algebraicRules):
                raise ValueError('System appears under-determined. '
                                 'Not enough  algebraic rules.') 
            self._makeDiffEqRHS()
            for func in self._dynamic_structure_funcs:
                exec 'self.%s_functionBody = self._make_%s()' % (func, func)
                _exec_dynamic_func(self, func, self.namespace)
            self._last_structure = copy.deepcopy(curr_structure)

        if self.events != getattr(self, '_last_events', None)\
           or structure_changed:
            logger.debug('Network %s: compiling events.' % self.id)
            for func in self._dynamic_event_funcs:
                exec 'self.%s_functionBody = self._make_%s()' % (func, func)
                _exec_dynamic_func(self, func, self.namespace)
            self._last_events = copy.deepcopy(self.events)

        # after compile, check that there are not more algebraic variables than
        # algebraic rules.
        if len(self.algebraicVars) > len(self.algebraicRules):
            raise ValueError('System appears under-determined. Not enough  algebraic rules.') 
        


    def _get_structure(self):
        """
        Return a tuple representing the structure of the Network.

        The tuple contains the functionDefinitions, reactions, and rules for
        the Network. It also contains information on the constancy and
        optimizability of the variables.
        """
        var_struct = {}
        structure = (self.functionDefinitions, self.reactions, 
                     self.assignmentRules, self.rateRules, var_struct,
                     self.algebraicRules)
        for id, var in self.variables.items():
            # If a constant variable is set equal to a function of other
            #  variables, we should include that function, otherwise
            #  our sensitivities will be wrong.
            if isinstance(self.get_var_ic(id), str):
                var_struct[id] = (var.is_constant, var.is_optimizable,
                                  self.get_var_ic(id))
            else:
                var_struct[id] = (var.is_constant, var.is_optimizable)

        return structure

    def takeDerivative(self, input, wrt):
        """
        Take the derivative of a math expression wrt a given variable id.

        Does the chain rule through assigned variables.
        """
        output = ExprManip.diff_expr(input, wrt)

        # What other assigned variables does input depend on?
        assigned_used = ExprManip.extract_vars(input)
        assigned_used.difference_update(sets.Set([wrt]))
        assigned_used.intersection_update(sets.Set(self.assignedVars.keys()))
        # Do the chain rule for those variables
        for id in assigned_used:
            rule = self.assignmentRules.getByKey(id)
            d2 = self.takeDerivative(rule, wrt)
            if d2 != '0':
                d = ExprManip.diff_expr(input, id)
                output += ' + (%s) *(%s)' % (d, d2)

        # What other constant variables does input depend on?
        constant_used = ExprManip.extract_vars(input)
        constant_used.difference_update(sets.Set([wrt]))
        constant_used.intersection_update(sets.Set(self.constantVars.keys()))
        # Do the chain rule for those variables
        for id in constant_used:
            ic = self.get_var_ic(id)
            if isinstance(ic, str):
                d2 = self.takeDerivative(ic, wrt)
                if d2 != '0':
                    d = ExprManip.diff_expr(input, id)
                    output += ' + (%s) *(%s)' % (d, d2)

        return ExprManip.simplify_expr(output)

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
        # deepcopy automatically does a deepcopy of whatever we return
        #  here, so we only need to do a shallow copy and remove functions 
        odict = copy.copy(self.__dict__)

        for func in self._dynamic_structure_funcs + self._dynamic_event_funcs:
            odict[func] = None
        odict['namespace'] = None
        # Let's not pickle these...
        odict['trajectory'] = None
        odict['ddv_dpTrajectory'] = None

        return odict

    def __setstate__(self, newdict):
        self.__dict__.update(newdict)
        self.namespace = copy.copy(self._common_namespace)

        # Recreate our namespace
        for func_id, func_str in self._func_strs:
            self.namespace[func_id] = eval(func_str, self.namespace, {})

        self._makeCrossReferences()

        # exec all our functions
        for func in self._dynamic_funcs:
            try:
                _exec_dynamic_func(self, func, self.namespace)
            except AttributeError:
                pass


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

    def get_eqn_structure(self):
        out = {}
        out['odes'] = dict(self.diff_eq_rhs.items())
        out['functions'] = {}
        for func_id, func_def in self.functionDefinitions.items():
            vars = ', '.join(func_def.variables)
            out['functions']['%s(%s)' % (func_id, vars)] = func_def.math
        out['parameters'] = dict([(id, var.value) for (id, var) 
                                  in self.constantVars.items()])
        out['assignments'] = dict(self.assignmentRules.items())
        out['events'] = dict([(event.trigger, event.event_assignments)
                              for event in self.events])

        return out


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

    def addEvent(self, id, trigger, eventAssignments, delay=0, name=''):
        self.add_event(id = id, trigger = trigger, 
                       event_assignments = eventAssignments, name = name,
                       delay = delay)

    addFunctionDefinition = add_func_def
    addAssignmentRule = add_assignment_rule

def _exec_dynamic_func(obj, func, in_namespace={}, bind=True):
    """
    Create the executable function corresponding to func's functionBody.
    """
    function_body = getattr(obj, '%s_functionBody' % func)
    # This exec gives the function access to everything defined in in_namespace
    #  and inserts the result into the locals namespace
    exec function_body in in_namespace, locals()
    # The call to types.MethodType ensures that we can call the function
    #  as obj.f(...) and get the implicit 'self' argument.
    # locals()[func] just gets the actual function object the exec created.
    #  Note that this this does depend on the _functionBody using a def
    #  with the proper name.
    setattr(obj, func, 
            types.MethodType(locals()[func], obj, obj.__class__))

    if HAVE_PSYCO and bind:
        psyco.bind(getattr(obj, func))

