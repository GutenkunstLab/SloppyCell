# This makes integer like one would expect mathematically, e.g. 1/2 = .5
#  rather than 0. This is obviously safer when loading other folks' models.
#  It does, however, cost about 10% performance on PC12.
# This is expected to become standard behavior around python 3.0
from __future__ import division

import copy
import sets
import types
import time
import os
import sys
import operator

import logging
logger = logging.getLogger('ReactionNetworks.Network_mod')

import scipy
import math

import SloppyCell
import SloppyCell.Utility as Utility
import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList

import SloppyCell.ExprManip as ExprManip
# We load a dictionary of previously-taken derivatives for efficiency
ExprManip.load_derivs(os.path.join(SloppyCell._TEMP_DIR, 'diff.pickle'))
# This will save the diffs dictionary upon exit from the python interpreter,
#  but only on the master node.
if SloppyCell.my_rank == 0:
    import atexit
    atexit.register(ExprManip.save_derivs, os.path.join(SloppyCell._TEMP_DIR, 
                                                        'diff.pickle'))

import Reactions
from Components import *
import Dynamics
import Trajectory_mod

_double_max_ = scipy.finfo(scipy.float_).max
_double_tiny_ = scipy.finfo(scipy.float_).tiny

class Network:
    # These are the methods that should be called to generate all the dynamic
    # functions. Each method should store the resulting function body in
    # self._dynamic_funcs_python or self._dynamic_funcs_c. C functions should
    # also have prototypes stored in self._prototypes_c.
    _dynamic_structure_methods = ['_make_res_function',
                                  '_make_alg_deriv_func',
                                  '_make_alg_res_func',
                                  '_make_dres_dc_function',
                                  '_make_dres_dcdot_function',
                                  '_make_ddaskr_jac',
                                  '_make_dres_dsinglep',
                                  '_make_sens_rhs',
                                  '_make_log_funcs',
                                  '_make_integrate_stochastic_tidbit'
                                  ]
    _dynamic_event_methods = ['_make_root_func']
    
    # These are-predefined functions we want in our working namespace
    _common_namespace = {'log': scipy.log,
                         'log10': scipy.log10,
                         'exp': math.exp,
                         'cos': math.cos,
                         'sin': math.sin,
                         'tan': math.tan,
                         'acos': math.acos,
                         'asin': math.asin,
                         'atan': math.atan,
                         'cosh': math.cosh,
                         'sinh': math.sinh,
                         'tanh': math.tanh,
                         # These don't have C support, so I've dropped them for
                         # now.
                         #'arccosh': scipy.arccosh,
                         #'arcsinh': scipy.arcsinh,
                         #'arctanh': scipy.arctanh,
                         'pow': math.pow,
                         'sqrt': math.sqrt,
                         'exponentiale': math.e,
                         'pi': math.pi,
                         'scipy': scipy,
                         'operator': operator,
                         'min': scipy.minimum,
                         'max': scipy.maximum
                         }
    # These are functions we need to create but that should be used commonly
    _standard_func_defs = [('root', ['n', 'x'],  'x**(1./n)'),
                           ('cot', ['x'],  '1./tan(x)'),
                           ('arccot', ['x'],  'atan(1/x)'),
                           ('coth', ['x'],  '1./tanh(x)'),
                           #('arccoth', ['x'],  'arctanh(1./x)'),
                           ('csc', ['x'],  '1./sin(x)'),
                           ('arccsc', ['x'],  'asin(1./x)'),
                           ('csch', ['x'],  '1./sinh(x)'),
                           #('arccsch', ['x'],  'arcsinh(1./x)'),
                           ('sec', ['x'],  '1./cos(x)'),
                           ('arcsec', ['x'],  'acos(1./x)'),
                           ('sech', ['x'],  '1./cosh(x)'),
                           #('arcsech', ['x'],  'arccosh(1./x)'),
                           ]
    # These define all the functions needed to deal with triggers.
    _logical_comp_func_defs = [('gt', ['x', 'y'],  'x > y'),
                               ('geq', ['x', 'y'],  'x >= y'),
                               ('lt', ['x', 'y'],  'x < y'),
                               ('leq', ['x', 'y'],  'x <= y'),
                               ('eq', ['x', 'y'],  'y == x'),
                               ('not_func', ['x'], 'not x'),
                               ('and_func', '*', 'x and y'),
                               ('or_func', '*', 'x or y'),
                               ]
    # Note that and_func and or_func are special-cased below, because they can
    # take a variable number of arguments in SBML. This has to be done because
    # the comp_func_defs list is used two ways. 1) It is used to substitute
    # within expressions derived from libsbml. 2) It is used to generate lambda
    # expressions for the functions.
                           
    # Add our common function definitions to a func_strs list.
    _common_func_strs = KeyedList()
    for id, vars, math in _standard_func_defs:
        var_str = ','.join(vars)
        func = 'lambda %s: %s' % (var_str, math)
        _common_func_strs.set(id, func)
        # These are all the partial derivatives
        for ii, wrt in enumerate(vars):
            deriv_id = '%s_%i' % (id, ii)
            diff_math = ExprManip.diff_expr(math, wrt)
            func = 'lambda %s: %s' % (var_str, diff_math)
            _common_func_strs.set(deriv_id, func)

    # Now add the C versions.
    _common_func_strs_c = KeyedList()
    _common_prototypes_c = KeyedList()
    for id, vars, math in _standard_func_defs:
        vars_c = ['double %s' % var for var in vars]
        var_str_c = ','.join(vars_c)
        c_body = []
        c_body.append('double %s(%s){' % (id, var_str_c))
        c_body.append('return %s;' % ExprManip.make_c_compatible(math))
        c_body.append('}')
        c_body = os.linesep.join(c_body)
        _common_prototypes_c.set(id, 'double %s(%s);' % (id, var_str_c))
        _common_func_strs_c.set(id, c_body)
        for ii, wrt in enumerate(vars):
            deriv_id = '%s_%i' % (id, ii)
            diff_math = ExprManip.diff_expr(math, wrt)
            c_body = []
            c_body.append('double %s(%s){' % (deriv_id, var_str_c))
            c_body.append('return %s;' % ExprManip.make_c_compatible(diff_math))
            c_body.append('}')
            c_body = os.linesep.join(c_body)
            _common_prototypes_c.set(deriv_id,'double %s(%s);' % (deriv_id, 
                                                                  var_str_c))
            _common_func_strs_c.set(deriv_id, c_body)

    # Also do the logical functions. These don't have derivatives.
    for id, vars, math in _logical_comp_func_defs:
        # and_func and or_func are special-cased, because the SBML versions can
        # take multiple arguments.
        if id == 'and_func':
            func = 'lambda *x: reduce(operator.and_, x)'
        elif id == 'or_func':
            func = 'lambda *x: reduce(operator.or_, x)'
        else:
            var_str = ','.join(vars)
            func = 'lambda %s: %s' % (var_str, math)
        _common_func_strs.set(id, func)

    # These are the strings we eval to create the extra functions our
    # Network defines. We'll evaluate them all to start with, to get the
    # common ones in there.
    for func_id, func_str in _common_func_strs.items():
        _common_namespace[func_id] = eval(func_str, _common_namespace, {})

    def __init__(self, id, name=''):
        self.id, self.name = id, name

        self.functionDefinitions = KeyedList()
        self.reactions = KeyedList()
        self.assignmentRules = KeyedList()
        self.rateRules = KeyedList()
        self.algebraicRules = KeyedList()
        self.events = KeyedList()
        self.constraints = KeyedList()        

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

        # All expressions are evaluated, and functions exec'd in self.namespace.
        # (A dictionary mapping names to objects they represent)
        self.namespace = copy.copy(self._common_namespace)
        self._func_strs = copy.copy(self._common_func_strs)

        # Option to integrate with log concentrations (to avoid negative
        # concentrations)
        self.integrateWithLogs = False

        self.len_root_func = 0

        # These dictionaries are for storing the function bodies of our dynamic
        # functions.
        self._dynamic_funcs_python = KeyedList()
        self._prototypes_c = self._common_prototypes_c.copy()
        self._dynamic_funcs_c = KeyedList()
        self._func_defs_c = KeyedList()

        self.deriv_funcs_enabled = True
        self.compiled = False

        self._manualCrossReferences_flag = False

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
        if not self._manualCrossReferences_flag:
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
                    is_optimizable=False, uniprot_ids=None):
        """
        Add a species to the Network.
        """
        species = Species(id, compartment, initial_conc, name, typical_value,
                          is_boundary_condition, is_constant, is_optimizable,
                          uniprot_ids)
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

    def add_event(self, id, trigger, event_assignments={}, delay=0, name='',
                  buffer=0):
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
        event = Event(id, trigger, event_assignments, delay, name,
                      buffer)
        self._checkIdUniqueness(event.id)
        self.events.set(event.id, event)

    def add_constraint(self, id, trigger, message=None, name=''):
        """
        Add a constraint to the Network.

        id - id for this constraint
        trigger - We treat constraints as events that correspond to an
            invalid solution whenever the trigger is True.
            Example: To have an invalid solution when species A is > 5.0:
                       trigger = 'lt(A, 5.0)'
        name - A more detailed name for the constraint, not restricted to the id
            format
        """
        constraint = ConstraintEvent(id, trigger, message, name)
        self._checkIdUniqueness(constraint.id)
        self.constraints.set(constraint.id, constraint)
        

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
        # Also do the evaluation.
        var_str = ','.join(variables)
        func_str = 'lambda %s: %s' % (var_str, math)
        self._func_strs.set(id, func_str)
        self.namespace[id] = eval(func_str, self.namespace)

    def make_func_defs(self):
        self.namespace = copy.copy(self._common_namespace)
        for id, func in self.functionDefinitions.items():
            variables = func.variables
            math = func.math

            var_str = ','.join(variables)
            func_str = 'lambda %s: %s' % (var_str, math)
            self._func_strs.set(id, func_str)
            logger.debug('Adding function to namespace: %s, %s' % (id, func_str))
            self.namespace[id] = eval(func_str, self.namespace)

            # Add the function and its partial derivatives to the C code.
            c_vars = ['double %s' % var for var in variables]
            c_var_str = ','.join(c_vars)
            self._prototypes_c.set(id,'double %s(%s);' % (id, c_var_str))
            c_math = ExprManip.make_c_compatible(math)
            c_body = 'double %s(%s){return %s;}' % (id, c_var_str, c_math)
            self._func_defs_c.set(id, c_body)

            for ii, wrt in enumerate(variables):
                # Python derivatives
                diff_id = '%s_%i' % (id, ii)
                diff_math = ExprManip.diff_expr(math, wrt)
                func_str = 'lambda %s: %s' % (var_str, diff_math)
                self._func_strs.set(diff_id, func_str)
                self.namespace[diff_id] = eval(func_str, self.namespace)

                # C derivatives
                self._prototypes_c.set(diff_id, 
                                       'double %s(%s);' % (diff_id,c_var_str))
                c_math = ExprManip.make_c_compatible(diff_math)
                c_body = 'double %s(%s){return %s;}' % (diff_id, c_var_str, 
                                                        c_math)
                self._func_defs_c.set(diff_id, c_body)

    def addReaction(self, id, *args, **kwargs):
        # Reactions can be added by (1) passing in a string representing
        #  kinetic law, or (2) passing in a class already specifying the 
        #  kinetic law.
        # XXX: I'm a little unhappy with this because option (2) breaks the
        #      pattern that the first argument is the id
        if isinstance(id, str):
            rxn = apply(Reactions.Reaction, (id,) + args, kwargs)
        else:
            rxn = apply(id, args, kwargs)

        self._checkIdUniqueness(rxn.id)
        self.reactions.set(rxn.id, rxn)

    def add_assignment_rule(self, var_id, rhs, index=None):
        """
        Add an assignment rule to the Network.

        A rate rules species that <var_id> = rhs.

        index: Optionally specify which index in the list of rules should be
               used for the new rule. This is important because assignment
               rules must be evaluated in order. The default is to add to the
               end of the list.
        """
        self.set_var_constant(var_id, False)
        if index is None:
            self.assignmentRules.set(var_id, rhs)
        else:
            self.assignmentRules.insert_item(index, var_id, rhs)
        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()
        # Put this assignment into effect...
        # The time=0 is somewhat arbitrary, reflecting the fact that we usually
        #  start integrations from time=0.
        # updateAssignedVars needs cross-references made, so we delay evaluation
        # if we're doing _manualCrossReferences
        if not self._manualCrossReferences_flag:
            self.updateAssignedVars(time=0)

    def add_rate_rule(self, var_id, rhs):
        """
        Add a rate rule to the Network.

        A rate rules species that d <var_id>/dt = rhs.
        """
        self.set_var_constant(var_id, False)
        self.rateRules.set(var_id, rhs)
        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()

    def add_algebraic_rule(self, rhs):
        """
        Add an algebraic rule to the Network.

        An algebraic rule specifies that 0 = rhs.
        """
        self.algebraicRules.set(rhs, rhs)
        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()

    def remove_component(self, id):
        """
        Remove the component with the given id from the Network.

        Components that can be removed are variables, reactions, events,
        function definitions, assignment rules, rate rules, and constraints.
        """
        complists = [self.variables, self.reactions, self.functionDefinitions,
                     self.events, self.constraints, self.assignmentRules, self.rateRules,
                     self.algebraicRules]
        for complist in complists:
            # If the id is in a list and has a non-empty name
            if complist.has_key(id):
                complist.remove_by_key(id)

        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()

    def _checkIdUniqueness(self, id):
        """
        Check whether a given id is already in use by this Network.
        """
        if id == 'time':
            logger.warn("Specifying 'time' as a variable is dangerous! Are you "
                        "sure you know what you're doing?")
        elif id == 'default':
            logger.warn("'default' is a reserved keyword in C. This will cause "
                        "problems using the C-based integrator.")
        elif id[0].isdigit():
            raise ValueError("The id %s is invalid. ids must not start with a "
                             "number." % id)
        if id in self.variables.keys()\
           or id in self.reactions.keys()\
           or id in self.functionDefinitions.keys()\
           or id in self.events.keys()\
           or id in self.constraints.keys()\
           or id == self.id:
            raise ValueError, ('The id %s is already in use!' % id)

    def set_id(self, id):
        """
        Set the id of this Network
        """
        if id != self.id:
            self._checkIdUniqueness(id)
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
        self.name = name 

    def get_name(self):
        """
        Get the name of this Network.
        """
        return self.name

    def set_deterministic(self):
        """
        Disables the stochastic simulation of network dynamics
        """
        if hasattr(self, 'stochastic'): delattr(self, 'stochastic')
        
    def set_stochastic(self, seed=None, fill_dt=None, rmsd=None):
        """
        set_stochastic enables the stochastic simulation of the network
        dynamics (instead of deterministic integration). This will discretize
        the dynamic variables in the model (making them integers). This
        simulation has only limited support for events*, and does not
        implement algebraic rules** or rates rules**. We use a kinetic MC
        algorithm (aka Gillespie algorithm). The complementary function
        'set_deterministic' may be used to disable stochastic simulations of
        network dynamics. Please note the results returned for exactly the
        times asked are interpolations between the surrounding points, thus
        concentrations may be floats instead of integers.

        * Since there is no exact way to implement events in a stochastic
          simulation, only a very primative method is present. This algorithm
          currently only supports time triggered events (more complex triggers
          may be added in the future) and fires them as soon after that time
          is passed as possible. Please note that it is the users
          responsibility to ensure the events are firing in a sensible
          way. Non-sensible output will arise when the time steps are too
          large and events are not fired near their target trigger times. The
          fired times are logged in the logger at the debug level.
        
        ** Possible future implementation, although parsing rate rules have the
           possible issue of creating negative reaction rates (which have no
           physical meaning). Add reactions instead of rate rules for stochastic
           simulations.
        
        Inputs:
        (seed=None, fill_dt=None, rmsd=None)
        seed -- RNG seed used in the simulation. If none, generated from the
                system time and os.urandom (if available).
        fill_dt -- Fill in the trajectory at the requested time interval. This
                   is the simplest way to fill in a trajectory.
        rmsd -- Maximum root mean squared distance in the dynamic variable
                space before taking a data point (a better trajectory filler).

        Outputs: None
        """
        if seed==None: 
            def rehash(key): # From numpy's mtrand module
                key += ~(key << 15)
                key ^=  (key >> 10)
                key +=  (key << 3)
                key ^=  (key >> 6)
                key += ~(key << 11)
                key ^=  (key >> 16)
                return key

            t = time.time()
            ts = int(t)
            tms = int((t-float(ts))*1000000.0)
            seed = rehash(ts)^rehash(tms)

            try:
                import struct
                seed ^= rehash(
                    struct.unpack('i', os.urandom(struct.calcsize('i')))[0])
            except: pass

            seed = int(seed % sys.maxint)
            
        if rmsd==None:
            rmsd = _double_max_
            
        self.stochastic = {'seed':seed, 'reseed':True,
                           'fill_dt':fill_dt, 'rmsd':rmsd}

        if hasattr(self, 'periodic'):
            logger.warn('Cannot find periodic solutions with a stochastic solution, disabling periodic solution finding')
            delattr(self, 'periodic')
            self.remove_component('_tau')
            
    def set_periodic(self, period, xtol=0.001, maxfun=100, phase=None,
                     minVel=None, level=2, log=False, rel=False):
        """
        set_periodic ensures that a period solution has satisfactorily
        approached a limit cycle before returning the solution. Please note that
        this has not been implemented for sensitivity calculations.

        Inputs:
        (period, xtol=0.001, maxfun=100, phase=None, minVel=None, level=2,
        log=False, rel=False)
        period -- initial guess of the period (required)
        xtol -- tolerance required for convergence (all chemicals). If set too
                high, numerical (in)precision may cause no stable limit cycle.
        maxfun -- maximum number of iterations (one trajectory found in each)
        phase -- set to a ('name', value) tuple to set the t=0 solution to a
                 particular state in systems which no forcing. Value may be
                 'max','min', or a float concentration (careful with the later).
        minVel -- used to detect fixed points:
                  when the (max IC velocity)**2 < minVel the integration
                  is halted, the stability set to True, and the period set to
                  0.0 (indicating nonperiodic, fixed point soln).
        level -- correlation level (>=1) to use to find the period. For example,
                 if level=2, the difference between x(0)-x(1*tau) and
                 x(0)-x(2*tau) is used.
        log -- Set to True to use logarithms of state variables, constraining
               the solution to non-zero values
        rel -- Set to True to use relative differences (x1-x0)/min(x0,x1) in
               the limit cycle period search.
               
        Outputs: (stored in net.periodic dictionary)
        tol -- tolerance of the solution (rmsd of ic1 to ic0 for each level)
        period -- actual period of oscillation
        stable -- whether a stable solution was found
        feval -- number of iterations
        stableIC -- stable initial condition

        Additionally, the period is accessible during the calculation and
        afterwards via the '_tau' variable (useful for assignments which
        rely on the period).
        """
        self.periodic = {'xtol': xtol, 'maxfun':maxfun, 'phase':phase,
                         'level':int(level), 'log':log, 'rel':rel,
                         'feval': 0, 'period': period, 'minVel': minVel,
                         'stable': False, 'stableIC': {}}
        if hasattr(self, 'stochastic'):
            logger.warn('Cannot find periodic solutions with a stochastic solution, switching to dynamic solution')
            self.set_deterministic()
        self.add_parameter('_tau', period, is_constant=True,
                           is_optimizable=False)
        
    def _iter_limit_cycle(self, params, varNames, s0):
        """
        Internal function used to integrate a trajectory and find the point
        nearest the initial point (defined as the root sum squared distance)
        in the period range of [0.475*tau, 1.5*tau]. The point found and period
        (tau) is returned.

        Inputs:
        (params, varNames, s0)
        params -- Network parameters
        varNames -- Names corresponding to the values in ln_x0
        s0 -- Initial state (period and conditions) where s0[0]=period and
              s0[1:]=initial conditions

        Outputs:
        (s1)
        s1 -- New state found (period and conditions defined same as input)
        """
        if self.periodic['rel']:
            rmsd = lambda _s1, _s2: \
                   scipy.sqrt(scipy.average(((_s1-_s2) / \
                                             scipy.minimum(_s1,_s2))**2))
        elif self.periodic['log']:
            rmsd = lambda _s1, _s2: \
                   scipy.sqrt(scipy.average((scipy.log(_s1/_s2)**2)))
        else:
            rmsd = lambda _s1, _s2: scipy.sqrt(scipy.average((_s1-_s2)**2))

        # Set the initial condition
        s0 = scipy.array(s0)
        tau0, x0 = s0[0], s0[1:]
        for name, value in zip(varNames, x0):
            self.set_var_ic(name, value)

        # Set the period search interval, restrictions:
        # 1) first guess by fminbound is at tau0
        # 2) Lower bound < 0.5*tau0 (to possibly detect double periods)
        # 3) Upper bound < 2.0*tau0 (to prevent double periods)
        a = 0.475
        b = a + (1.0-a)/(0.5*(3.0-scipy.sqrt(5.0))) # ~1.85
        
        # Find the trajectory
        max_time = b * tau0 * float(self.periodic['level'])
        self.trajectory = self.integrate([0.,max_time], params, addTimes = True)
        self.trajectory.build_interpolated_traj()
        ev_inter_traj = self.trajectory.evaluate_interpolated_traj
        
        # Find the period between 1/2*tau0 and 3/2*tau0
        # Without the squeeze, this failed on some versions of numpy when 
        # t = array([0.2])
        s = lambda t: scipy.array([scipy.squeeze(t)] + \
                                  [ev_inter_traj(name, t) for name in varNames])
        
        def __iter_cost(t):
            t = scipy.squeeze(abs(t))
            res = scipy.array([t], scipy.float_)
            for lvl in range(1, self.periodic['level']+1):
                res = scipy.concatenate((res, s(t*float(lvl))[1:]))
            res0 = scipy.concatenate(([tau0], list(x0)*self.periodic['level']))
            return rmsd(res, res0)

        tau1 = scipy.optimize.fminbound(__iter_cost, a * tau0, b * tau0)
        tau1 = abs(tau1)
        s1 = s(tau1)

        # Update the periodic dictionary
        self.periodic['feval'] += 1
        self.periodic['period'] = tau1
        self.set_var_val('_tau', tau1)
        self.periodic['tol'] = __iter_cost(tau1)
        if self.periodic['tol'] < self.periodic['xtol']:
            self.periodic['stable']=True

        logger.debug('Limit Cycle Iteration #%i: period = %f, tol = %f, %s'% \
                     (self.periodic['feval'], self.periodic['period'],
                      self.periodic['tol'],
                      '%sstable'%('un'*(not self.periodic['stable']))))
        
        return s1

    def _eliminate_slowest_mode(self, s0, s1, s2):
        """
        Internal function which finds and eliminates the slowest decaying
        mode in the approach to a limit cycle:
        si = [tau_i, x0_0, .., x0_i]
        F0 = s_star + vector_1 + ...
        F1 = s_star + lambda_1 * vector_1 + ...
        F2 = s_star + lambda_1^2 * vector_1 + ...

        lambda_1 = (F2-F1)*(F1-F0)/|F1-F0|^2
        vector_1 = (F2-F1)/(lambda_1 * (lambda_1 - 1.0))

        Inputs: (s0, s1, s2)
        Three vectors representing three consecutive steps towards a limit cycle

        Ouptuts: (s_star)
        s_star = s2 - lambda_1^2 * vector_1
        """
        x0, x1, x2 = s0[1:], s1[1:], s2[1:]
        lambda_1  = sum( (x2-x1) * (x1-x0) ) / sum( (x1-x0)**2 )
        vector_1 = (x2-x1)/(lambda_1 * (lambda_1-1.))
        s = scipy.array(s2) # make a copy of s2
        s[1:] -= lambda_1*lambda_1*vector_1
        return s

    def _find_limit_cycle(self, params):
        """
        Internal function to find a stable limit cycle (if one exists), given a
        parameter set, an estimated period, and a set of initial conditions.
        This function calculates two iterations (which are hopefully approaching
        a limit cycle), and use this to estimate a point on the limit cycle - an
        estimation which is kept if it is closer to the limit cycle.

        Inputs:
        (params,)
        params -- model parameters

        Outputs: None
        """
        self.periodic['stable'] = False
        self.periodic['feval'] = 0
        
        # We only want dynamic variables whose value is not assigned
        varNames = [name for name in self.dynamicVars.keys()
                   if name in self.species.keys()]

        # The initial period and state
        tau0 = self.periodic['period']
        x0 = [self.periodic['stableIC'].get(name, self.get_var_ic(name))
              for name in varNames]
        s0 = scipy.array([tau0]+x0)

        while not self.periodic['stable'] and \
                  self.periodic['feval'] < self.periodic['maxfun']:

            # Perform the first two searched
            s1 = self._iter_limit_cycle(params, varNames, s0)
            if self.periodic['stable']:
                s0 = s1
                break

            s2 = self._iter_limit_cycle(params, varNames, s1)
            tol2 = self.periodic['tol']
            if self.periodic['stable']:
                s0 = s2
                break

            if self.periodic['log']:
                ls0, ls1, ls2 = scipy.log(s0), scipy.log(s1), scipy.log(s2)
                s_star = scipy.exp(self._eliminate_slowest_mode(ls0, ls1, ls2))
            else:
                s_star = self._eliminate_slowest_mode(s0, s1, s2)
            
            # Predict a point on the orbit, and choose the best point to
            # use as the next round (either the last search #2 or the
            # prediction)
            try:
                s3 = self._iter_limit_cycle(params, varNames, s_star)

                # Did we find a better solution?
                if self.periodic['tol'] < tol2: s0 = s3
                else: s0 = s2  # Prediction (s_star) was worse
            except:
                s0 = s2
                
            # Found NANs, cut losses and run
            if scipy.any(scipy.isnan(s0)):
                self.periodic['stable'] = False
                break

            # Test for a fixed point, stopping if one is found
            if self.periodic['minVel'] is not None:
                for name, val in zip(varNames, s0[1:]):
                    self.set_var_ic(name, val)

                maxVel = max(
                    [iV**2 for name, iV in zip(self.dynamicVars.keys(),
                                               self.get_initial_velocities())
                     if name in varNames])

                if maxVel < self.periodic['minVel']:
                    self.periodic['period']=0.0
                    self.set_var_val('_tau', 0.0)
                    self.periodic['stable']=True
                    break


        if self.periodic['stable']:
            # Record the initial conditions (if stable)
            for name, value in zip(varNames, s0[1:]):
                self.periodic['stableIC'][name] = value
                self.set_var_ic(name, value)

            # Set the phase of the cycle (or else return, we're done)
            if self.periodic['phase'] is not None:
                # Define the cost
                pN, pV = self.periodic['phase']
                _eit = self.trajectory.evaluate_interpolated_traj

                if pV=='min': cost = lambda t: _eit(pN, t)
                elif pV=='max': cost = lambda t: 1.0 / _eit(pN, t)
                else: cost = lambda t: (pV - _eit(pN, t))**2

                # Find the best-fitting point
                tau = self.periodic['period']
                lowerB = float(self.periodic['level'] - 1) * tau
                upperB = float(self.periodic['level']) * tau
                t0 = scipy.optimize.fminbound(cost, lowerB, upperB)

                # Set the initial point, fall through and integrate
                for name in varNames:
                    self.periodic['stableIC'][name] = _eit(name, t0)
                    self.set_var_ic(name, _eit(name, t0))


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
            elif var.endswith('_maximum') or var.endswith('_minimum'):
                t1,t2 = times
                if t1 is not None:
                    t.add(t1)
                if t2 is not None:
                    t.add(t2)
            elif self.variables.has_key(var):
                t.union_update(sets.Set(times))
            else:
                raise ValueError('Unknown variable %s requested from network %s'
                                 % (var, self.get_id()))

        # This takes care of the normal data points
        t = list(t)
        t.sort()

        if hasattr(self, 'periodic'):
            self._find_limit_cycle(params)
            if self.periodic['period'] > max(t):
                t.append(self.periodic['period'])
        elif hasattr(self, 'stochastic'):
            self.trajectory = self.integrateStochastic(t, params=params)
            return

        self.trajectory = self.integrate(t, params=params, 
                                         addTimes=ret_full_traj)
        
    def CalculateSensitivity(self, vars, params):
        t = sets.Set([0])

        for var,times in vars.items():
            if var.endswith('_maximum') or var.endswith('_minimum'):
                t1,t2 = times
                if t1 is not None:
                    t.add(t1)
                if t2 is not None:
                    t.add(t2)
            elif self.variables.has_key(var):
                t.union_update(sets.Set(times))
            else:
                raise ValueError('Unknown variable %s requested from network %s'
                                 % (var, self.get_id()))

        t = list(t)
        t.sort()

        self.ddv_dpTrajectory = self.integrateSensitivity(t,params=params, 
                                                          addTimes=True)
        self.trajectory = self.ddv_dpTrajectory

    def GetName(self):
        return self.get_id()

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
            if id.endswith('_maximum') or id.endswith('_minimum'):
                # This is an extremum search.
                # Which variable are we looking for?
                var = '_'.join(id.split('_')[:-1])

                # First, check though the trajectory, to see if the extremum
                # is in there.
                t = self.trajectory.get_times()
                minTime, maxTime = vars[id]
                minSearchTime, maxSearchTime = vars[id]
                if minSearchTime is None:
                    minSearchTime = t[0]
                if maxSearchTime is None:
                    maxSearchTime = t[-1]
                # Restrict timespan in which to search.
                min_t_ii = t.searchsorted(minSearchTime)
                max_t_ii = t.searchsorted(maxSearchTime)-1

                vartraj = self.trajectory.get_var_traj(var)
                if id.endswith('_maximum'):
                    var_val_ii = vartraj[min_t_ii:max_t_ii+1].argmax()\
                            + min_t_ii
                elif id.endswith('_minimum'):
                    var_val_ii = vartraj[min_t_ii:max_t_ii+1].argmin()\
                            + min_t_ii
                var_val = vartraj[var_val_ii]
                t_val = t[var_val_ii]

                # Now check through the events. If we're running under 
                # full_speed the trajectory can be very sparse. If the user has
                # set up events to tag potential extrema, this will use those.
                curr_vals = self.get_var_vals()

                isDynamicVar = self.dynamicVars.has_key(var)
                if isDynamicVar:
                    dynVarIndex = self.dynamicVars.index_by_key(var)

                for holder in self.trajectory.events_occurred:
                    possibilities = [(holder.time_fired, holder.y_fired),
                                     (holder.time_exec, holder.y_pre_exec),
                                     (holder.time_exec, holder.y_post_exec)]
                    # There are three possible event-related times. When the
                    # event fired, when it executed, and values pre- and post-
                    # execution.
                    for poss_ii, (time, dynVarVals) in enumerate(possibilities):
                        if isDynamicVar:
                            eval = dynVarVals[dynVarIndex]
                        elif hasattr(holder, 'all_var_vals')\
                                and len(holder.all_var_vals) >= poss_ii + 1:
                            eval = holder.all_var_vals[poss_ii].get(var)
                        else:
                            if not hasattr(holder, 'all_var_vals'):
                                holder.all_var_vals = []
                            self.updateVariablesFromDynamicVars(dynVarVals, 
                                                                time)
                            holder.all_var_vals.append(self.get_var_vals())
                            eval = self.get_var_val(var)

                        if id.endswith('_maximum') and eval > var_val\
                           and time >= minSearchTime and time <= maxSearchTime:
                            var_val = eval
                            t_val = time
                        elif id.endswith('_minimum') and eval < var_val\
                           and time >= minSearchTime and time <= maxSearchTime:
                            var_val = eval
                            t_val = time
                self.set_var_vals(curr_vals)

                result[id] = {(minTime, maxTime):(t_val,var_val)}
            elif self.variables.has_key(id):
                traj = self.trajectory.getVariableTrajectory(id)
                result[id] = dict(zip(times, traj))

        return result

    def GetSensitivityResult(self, vars):
       # note: returns all the timepoints we have, not just
       # the requested ones (which should be a subset of all
       # the timepoints)
       result = {}
       times = self.ddv_dpTrajectory.get_times()
       for id in vars.keys():
           result[id] = {}
           if id.endswith('_maximum') or id.endswith('_minimum'):
               # Sensitivity of extremum:
               # d/dp y(tmax) = dy/dp + dy/dt*dtmax/dp, but dy/dt should be zero
               # at a max, so it's just dy/dp...
               # This is an extremum search.
               # Which variable are we looking for?
               var = '_'.join(id.split('_')[:-1])

               # First, check though the trajectory, to see if the extremum
               # is in there.
               t = self.trajectory.get_times()
               minTime, maxTime = vars[id]
               minSearchTime, maxSearchTime = vars[id]
               if minSearchTime is None:
                   minSearchTime = t[0]
               if maxSearchTime is None:
                   maxSearchTime = t[-1]
               # Restrict timespan in which to search.
               min_t_ii = t.searchsorted(minSearchTime)
               max_t_ii = t.searchsorted(maxSearchTime)-1

               vartraj = self.trajectory.get_var_traj(var)
               if id.endswith('_maximum'):
                   var_val_ii = vartraj[min_t_ii:max_t_ii+1].argmax()\
                           + min_t_ii
               elif id.endswith('_minimum'):
                   var_val_ii = vartraj[min_t_ii:max_t_ii+1].argmin()\
                           + min_t_ii
               var_val = vartraj[var_val_ii]
               t_val = t[var_val_ii]
               var_val_sens = dict((optvar, self.ddv_dpTrajectory.get_var_val_index((var,optvar), var_val_ii)) for optvar in self.optimizableVars.keys())

               # Now check through the events. If we're running under 
               # full_speed the trajectory can be very sparse. If the user has
               # set up events to tag potential extrema, this will use those.
               curr_vals = self.get_var_vals()
               Ndynvars = len(self.dynamicVars)
               try:
                   var_index = self.dynamicVars.keys().index(var)
               except ValueError:
                   logger.warn("Can't check events for %s because %s is not a "
                               "dynamic variable. (It is probably assigned.)"
                               % (id, var))
                   result[id][(minTime, maxTime)] = var_val_sens
                   continue
               for holder in self.trajectory.events_occurred:
                   possibilities = [(holder.time_fired, holder.ysens_fired),
                                    (holder.time_exec, holder.ysens_pre_exec),
                                    (holder.time_exec, holder.ysens_post_exec)]
                   # There are three possible event-related times. When the
                   # event fired, when it executed, and values pre- and post-
                   # execution.
                   for time, ysens in possibilities:
                       self.updateVariablesFromDynamicVars(ysens[:Ndynvars], 
                                                           time)
                       eval = self.get_var_val(var)
                       if (time >= minSearchTime and time <= maxSearchTime)\
                          and ((id.endswith('_maximum') and eval > var_val)
                               or (id.endswith('_minimum') and eval < var_val)):
                           var_val = eval
                           t_val = time
                           var_val_sens = dict(zip(self.optimizableVars.keys(),
                                                   ysens[Ndynvars:][var_index::Ndynvars]))
               self.set_var_vals(curr_vals)

               result[id][(minTime, maxTime)] = var_val_sens
           for tIndex, t in enumerate(times):
               result[id][t] = {}
               for optparam in self.optimizableVars.keys():
                   result[id][t][optparam] = \
                       self.ddv_dpTrajectory.get_var_val_index((id,optparam), tIndex)
       return result

    def integrate(self, times, params=None, returnEvents=False, addTimes=True,
                  rtol = None):
        if self.add_tail_times or addTimes:
            times = scipy.concatenate((times, [1.05*times[-1]]))
        return Dynamics.integrate(self, times, params=params, 
                                  fill_traj=(self.add_int_times or addTimes))

    def integrateSensitivity(self, times, params = None,
                             returnEvents = False, addTimes = True,
                             rtol=None):
        if self.add_tail_times:
            times = scipy.concatenate((times, [1.05*times[-1]]))
        return Dynamics.integrate_sensitivity(self, times, params, rtol,
                                              fill_traj=self.add_int_times)

    def integrateStochastic(self, times, params=None):
        if self.stochastic['fill_dt'] is not None:
            times = sets.Set(times)
            times.union_update(scipy.arange(min(times), max(times),
                                            self.stochastic['fill_dt']))
            times = list(times)
            times.sort()

        if times[0]==0.: # Only reset if specifically requested (not by event)
            self.resetDynamicVariables()

        # Add in the event times (only if they don't extend the trajectoy!)
        times = sets.Set(times)
        times.union_update([_.triggeringTime for _ in self.events
                            if _.triggeringTime<max(times)])
        times = list(times)
        times.sort()
        
        if times[0]==0.:
            times = times[1:]
            
        if params is not None: self.update_optimizable_vars(params)

        self.compile()

        # Test that we have a working integrator
        body = self._dynamic_funcs_python.get('integrate_stochastic_tidbit')
        if 'pass' in body.split('\n')[-1]:
            err_body = '\n'.join(
                ['Integrate stochastic failed due to the problem(s):']+
                body.split('\n')[1:-1])
            raise RuntimeError(err_body)
        
        dv=scipy.array([self.get_var_val(_) for _ in self.dynamicVars.keys()])
        cv = self.constantVarValues

        trajectory = Trajectory_mod.Trajectory(self, holds_dt=0, const_vals=cv)
        trajectory.appendFromODEINT(scipy.array([0.0]), scipy.array([dv]))

        # VERY simple event handling...
        events_occurred = []
        pendingEvents = {}
        for event in self.events.values():
            pendingEvents.setdefault(event.triggeringTime, [])
            pendingEvents[event.triggeringTime].append(event)
        
        t, tInd = 0.0, 0
        while t < times[-1]:
            # Get the next chunk
            t, dv, tout, dvout = self.integrate_stochastic_tidbit(
                self.stochastic['seed'], self.stochastic['reseed'],
                t, dv, cv, self.stochastic['rmsd'], times[tInd])
            self.stochastic['reseed']=False

            # Fill in trajectory and execute any events that occurred
            while tInd<len(times) and t >= times[tInd]:
                # Extrapolate the dynamic vars, if not the stop_time requested
                if times[tInd]!=tout:
                    dt = t - trajectory.timepoints[-1]
                    dvout = trajectory.values[-1]+(dv-trajectory.values[-1])/\
                            dt*(times[tInd] - trajectory.timepoints[-1])

                trajectory.appendFromODEINT(scipy.array([times[tInd]]),
                                            scipy.array([dvout]))

                # Is this where an event fires?
                if times[tInd] in pendingEvents.keys():
                    events = pendingEvents.pop(times[tInd])
                    for event in events:
                        # Execute twice b/c we want to update the dvout at
                        # the fired time in addition to updating the dv for
                        # any future extrapolations. This is really a product
                        # of not being able to ensure we haven't overstepped
                        # our desired timepoint and crossed the next.
                        dvout = self.executeEvent(event, times[tInd], dvout, \
                                                  dvout, times[tInd])
                        dv = self.executeEvent(event, times[tInd], dvout, dv, t)
                        events_occurred.append(event)
                        logger.debug('Executed event %s in net %s at t=%s'%\
                                     (event.id, self.id, t))
                    trajectory.appendFromODEINT(scipy.array([times[tInd]]),
                                                scipy.array([dvout]))
                    
                tInd += 1

            if trajectory.timepoints[-1]<t:
                trajectory.appendFromODEINT(scipy.array([t]),
                                            scipy.array([dv]))

        return trajectory
    
    def _sub_for_piecewise(self, expr, time):
        """
        Runs through expr and replaces all piecewise expressions by the
        clause that is currently active.
        """
        if not expr.count('piecewise('):
            return expr
        funcs_used = ExprManip.extract_funcs(expr)
        for func, args in funcs_used:
            if func == 'piecewise':
                ast = ExprManip.strip_parse(expr)
                ast = self._sub_for_piecewise_ast(ast, time)
                return ExprManip.ast2str(ast)
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

    def evaluate_expr(self, expr, time=0, var_vals=None):
        """
        Evaluate the given expression using the current values of the network
        variables.
        """
        try:
            return float(expr)
        except ValueError:
            # We create a local_namespace to evaluate the expression in that
            #  maps variable ids to their current values
            expr = self._sub_for_piecewise(expr, time)
            if var_vals is None:
                vars_used = ExprManip.extract_vars(expr)
                var_vals = [(id, self.get_var_val(id)) for id in vars_used
                            if id != 'time']
                var_vals = dict(var_vals)
                var_vals['time'] = time
            local_namespace = var_vals
            local_namespace.update(self.namespace)
            # We strip whitespace, just for convenience
            return eval(expr.strip(), local_namespace, {})

    #
    # Methods to get and set object properties
    #
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
        return self.get_variable(id).typicalValue

    def get_var_typical_vals(self):
        """
        Return the variable typical values as a KeyedList.
        """
        return KeyedList([(id, self.get_var_typical_val(id)) for id in
                          self.variables.keys()])

    def set_var_ic(self, id, value, warn=True, update_constants=True):
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
            if update_constants:
                self.constantVarValues = [self.evaluate_expr(var.value) for var
                                          in self.constantVars.values()]
                self.constantVarValues = scipy.array(self.constantVarValues)

    def get_var_ic(self, id):
        """
        Return the initial condition for a variable
        """
        return self.get_variable(id).initialValue

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
        val = self.get_variable(id).value
        return self.evaluate_expr(val)

    def get_var_vals(self, ids=None):
        """
        Return the current variable values as a KeyedList

        ids -- List o variable ids to return values for. If None, all variables
               are used.
        """
        if ids is None:
            ids = self.variables.keys()
        return KeyedList([(id, self.get_var_val(id)) for id in ids])

    def get_initial_velocities(self):
        """
        Returns the vector field evaluated at the initial conditions
        """
        ics = [self.evaluate_expr(self.getInitialVariableValue(dvid))
               for dvid in self.dynamicVars.keys()]
        # the evaluate_expr is in case an initial condition is a parameter
        y0, initialv = Dynamics.find_ics(y=ics, yp=scipy.ones(len(ics)), 
                                         time=0, 
                                         var_types=self._dynamic_var_algebraic,
                                         rtol=[1e-6]*len(ics), 
                                         atol=[1e-6]*len(ics),
                                         constants=self.constantVarValues, 
                                         net=self)
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
                self.set_var_ic(id, params.get(id), update_constants=False)
        elif len(params) == len(self.optimizableVars):
            for ii, id in enumerate(self.optimizableVars.keys()):
                self.set_var_ic(id, params[ii], update_constants=False)
        else:
            raise ValueError('Passed in parameter set does not have the proper '
                             'length!')

        self.constantVarValues = [self.evaluate_expr(var.value) for var in
                                  self.constantVars.values()]
        self.constantVarValues = scipy.array(self.constantVarValues)

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
            if isinstance(var.initialValue, str):
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
        to_get = sets.Set(self.variables.keys())
        to_get.difference_update(sets.Set(self.assignedVars.keys()))
        var_vals = [(id, self.get_var_val(id)) for id in to_get]
        var_vals = dict(var_vals)
        var_vals['time'] = time
        for id, rhs in self.assignmentRules.items():
            assigned_val = self.evaluate_expr(rhs, time, var_vals=var_vals)
            self.assignedVars.getByKey(id).value = assigned_val
            var_vals[id] = assigned_val

    def set_var_optimizable(self, id, is_optimizable):
        self.get_variable(id).is_optimizable = is_optimizable
        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()

    def set_var_constant(self, id, is_constant):
        self.get_variable(id).is_constant = is_constant
        if not self._manualCrossReferences_flag:
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

    def _make_res_function(self):
        py_body = []
        py_body.append('def res_function(time, dynamicVars, yprime, '
                       'constants):')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yprime = scipy.asarray(yprime)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('residual = scipy.empty(%i, scipy.float_)'
                       % len(self.dynamicVars))
        py_body.append('')
        self._add_assignments_to_function_body(py_body)
        py_body.append('')

        c_body = []
        # The C version needs to take a number of arguments that the f2py
        # wrapper will hide. 
        # To be callable from Fortran it also needs 1) its name followed by an
        # underscore, and 2) to take all arguments as pointers.
        c_args = 'double *time_ptr, double *dynamicVars, double *yprime, '\
                'double *cj_ptr, double *residual, int *ires_ptr, '\
                'double *constants, int *ipar'
        c_body.append('void res_function_(%s){' % c_args)
        self._prototypes_c.set('res_function', 
                               'void res_function_(%s);' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True)
        c_body.append('')
                      
        # Make a list of algebraic rules for accessing later
        algebraicRuleList = self.algebraicRules.values()

        # We keep a list of terms in our residual function for use building the 
        #  various derivative functions.
        self._residual_terms = []
        # This list tells us whether a give dynamic variable is algebraic or not
        # It will be used in the integration.
        self._dynamic_var_algebraic = []

        # Loop over everything in the dynamicVars list
        for ii, (id, var) in enumerate(self.dynamicVars.items()):
            if self.algebraicVars.has_key(id):
                # It's an algebraic equation. Pop the first algebraic rule off
                # the list.
                rhs = algebraicRuleList.pop(0)
                py_body.append('# Residual function corresponding to an '
                               'algebraic  equation')
                py_body.append('residual[%i] = %s' % (ii, rhs))
                c_rhs = ExprManip.make_c_compatible(rhs)
                c_body.append('residual[%i] = %s;' % (ii, c_rhs))
                self._residual_terms.append((rhs, None))
                self._dynamic_var_algebraic.append(-1)
            else:
                rhs = self.diff_eq_rhs.getByKey(id)
                c_rhs = ExprManip.make_c_compatible(rhs)
                self._dynamic_var_algebraic.append(+1)
                py_body.append('# Residual function for %s' % id)
                if rhs != '0':
                    py_body.append('residual[%i] = %s - yprime.item(%i)' % 
                                (ii, rhs, ii))
                    c_body.append('residual[%i] = %s - yprime[%i];' % 
                                  (ii, c_rhs, ii))
                    self._residual_terms.append((rhs, id))
                else:
                    py_body.append('residual[%i] = -yprime.item(%i)' % (ii, ii))
                    c_body.append('residual[%i] = -yprime[%i];' % (ii, ii))
                    self._residual_terms.append((None, id))

        py_body.append('')
        py_body.append('return residual')
        py_body = '\n    '.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('res_function', py_body)
        self._dynamic_funcs_c.set('res_function', c_body)
                    
        return py_body
        
    def _make_root_func(self):
        len_root_func = 0

        py_body = []
        py_body.append('def root_func(time, dynamicVars, yprime, constants):')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yprime = scipy.asarray(yprime)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        # We don't know the length of root_devs yet, so for now we'll
        #  just insert a placeholder line.
        py_body.append('root_devs = scipy.empty(NEED_TO_FIX, scipy.float_)')
        py_body.append('')
        self._add_assignments_to_function_body(py_body, include_dts=self.deriv_funcs_enabled)
        py_body.append('')

        c_body = []
        c_args = 'int *neq_ptr, double *time_ptr, double *dynamicVars, '\
                'double *yprime, int *nrt_ptr, double *root_devs, '\
                'double *constants, int *ipar'
        c_body.append('void root_func_(%s){' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True, include_dts=self.deriv_funcs_enabled)
        self._prototypes_c.set('root_func', 
                               'void root_func_(%s);' % c_args)
        c_body.append('')

        self.event_clauses = []
        for ii, event in enumerate(self.events.values()+self.constraints.values()):
            trigger = event.trigger
            for func_name, func_vars, func_expr in self._logical_comp_func_defs:
                trigger = ExprManip.sub_for_func(trigger, func_name, func_vars,
                                                 func_expr)

            py_body.append('root_devs[%i] = (%s) - 0.5' % (len_root_func, 
                                                           trigger))
            c_trigger = ExprManip.make_c_compatible(trigger)
            c_body.append('root_devs[%i] = (%s) - 0.5;' % (len_root_func, 
                                                           c_trigger))
            self.event_clauses.append(trigger)
            len_root_func += 1

        for ii, event in enumerate(self.events.values()+self.constraints.values()):
            event.sub_clause_indices = []
            trigger = event.trigger
            for func_name, func_vars, func_expr in self._logical_comp_func_defs:
                trigger = ExprManip.sub_for_func(trigger, func_name, func_vars,
                                                 func_expr)
            comparisons = ExprManip.extract_comps(trigger)
            if len(comparisons) > 1:
                for comp in comparisons:
                    py_body.append('root_devs[%i] = (%s) - 0.5\n\t' 
                                   % (len_root_func, comp))
                    c_comp = ExprManip.make_c_compatible(comp)
                    c_body.append('root_devs[%i] = (%s) - 0.5;\n\t' 
                                   % (len_root_func, c_comp))
                    self.event_clauses.append(comp)
                    event.sub_clause_indices.append(len_root_func)
                    len_root_func += 1

        self.len_root_func = len_root_func
        # Now that we know how many roots we're looking for, go back and
        #  insert the proper size for root_devs
        py_body[5] = ('root_devs = scipy.empty(%i, scipy.float_)'
                      % len_root_func)
        py_body.append('')
        py_body.append('return root_devs')
        py_body = '\n\t'.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('root_func', py_body)
        self._dynamic_funcs_c.set('root_func', c_body)

        return py_body

    def _make_alg_deriv_func(self):
        # This function is used when we want to calculate consistent derivatives
        #  for algebraic variables. It returns the derivative wrt time of
        #  all the algebraic rules. Since the rules are <0=rule>, these should
        #  all be zero when we have consistent values for the variable derivs.
        # E.g. If we have an algebraic rule 0 = f(x,y), then taking the time
        #  derivative yields: 0 = df/dx*dx/dt + df/dy*dy/dt. This function
        #  returns df/dx*dx/dt + df/dy*dy/dt and fsolve twiddles with the
        #  dx/dt's that correspond to algebraic variables. It's a linear systems
        #  of equations, so fsolve is kind of overkill, but it hasn't been a
        #  slow point so far.
        py_body = []
        py_body.append('def alg_deriv_func(alg_yp, dynamicVars, yp, time, '
                       'constants):')
        py_body.append('alg_yp = scipy.asarray(alg_yp)')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yp = scipy.asarray(yp)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('alg_derivs_res = scipy.empty(%i, scipy.float_)'
                       % len(self.algebraicRules))
        py_body.append('')
        self._add_assignments_to_function_body(py_body)
        py_body.append('')

        c_body = []
        c_args = 'double *alg_yp, double *dynamicVars, double *yp, '\
                'double *time_ptr, double *constants, double *alg_derivs_res' 
        c_body.append('void alg_deriv_func_(%s){' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True)
        self._prototypes_c.set('alg_deriv_func', 
                               'void alg_deriv_func_(%s);' % c_args)
        c_body.append('')

        for rule_ii, rule in enumerate(self.algebraicRules.values()):
            # Each rhs is a sum of drule/dA * dA/dt + drule/dB * dB/dt + ...
            rhs_terms = []
            rhs_terms_c = []
            rule_vars = ExprManip.extract_vars(rule)
            deriv_wrt_time = self.takeDerivative(rule, 'time', rule_vars)
            rhs_terms.append(deriv_wrt_time)
            rhs_terms_c.append(deriv_wrt_time)
            for dyn_id in self.dynamicVars.keys():
                deriv = self.takeDerivative(rule, dyn_id, rule_vars)
                if deriv != '0':
                    if self.algebraicVars.has_key(dyn_id):
                        index = self.algebraicVars.index_by_key(dyn_id)
                        rhs_terms.append('(%s)*alg_yp.item(%i)' % (deriv,index))
                        rhs_terms_c.append('(%s)*alg_yp[%i]' % (deriv,index))
                    else:
                        index = self.dynamicVars.index_by_key(dyn_id)
                        rhs_terms.append('(%s)*yp.item(%i)' % (deriv, index))
                        rhs_terms_c.append('(%s)*yp[%i]' % (deriv, index))
            rhs = ' + '.join(rhs_terms)
            rhs_c = ' + '.join(rhs_terms_c)
            rhs_c = ExprManip.make_c_compatible(rhs_c)
            py_body.append('alg_derivs_res[%i] = %s' % (rule_ii, rhs))
            c_body.append('alg_derivs_res[%i] = %s;' % (rule_ii, rhs_c))

        py_body.append('')
        py_body.append('return alg_derivs_res')
        py_body = '\n\t'.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('alg_deriv_func', py_body)
        self._dynamic_funcs_c.set('alg_deriv_func', c_body)

        return py_body

    def _make_alg_res_func(self):
        py_body = []
        py_body.append('def alg_res_func(alg_vals, dynamicVars, '
                       'time, constants):')
        py_body.append('alg_vals = scipy.asarray(alg_vals)')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('residual = scipy.zeros(%i, scipy.float_)'
                       % len(self.algebraicVars))

        c_body = []
        c_args = 'double *alg_vals, double *dynamicVars, double *time_ptr, '\
                'double *constants, double *residual'
        self._prototypes_c.set('alg_res_func',
                               'void alg_res_func_(%s);' % c_args)
        c_body.append('void alg_res_func_(%s){' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')

        # Copy our algebraic variable guesses into the appropriate slots of
        # dynamicVars
        for ii, key in enumerate(self.algebraicVars.keys()):
            var_index = self.dynamicVars.index_by_key(key)
            py_body.append('dynamicVars[%i] = alg_vals[%i]'
                           % (var_index, ii))
            c_body.append('dynamicVars[%i] = alg_vals[%i];'
                          % (var_index, ii))

        # Now make all the assignments from dynamicVars to the variable ids
        py_body.append('')
        self._add_assignments_to_function_body(py_body)
        py_body.append('')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True)
        c_body.append('')

        # Now evaluate all the algebraic rules.
        for ii, rhs in enumerate(self.algebraicRules.values()):
            py_body.append('residual[%i] = %s' % (ii, rhs))
            c_rhs = ExprManip.make_c_compatible(rhs)
            c_body.append('residual[%i] = %s;' % (ii, c_rhs))

        py_body.append('')
        py_body.append('return residual')
        py_body = '\n    '.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('alg_res_func', py_body)
        self._dynamic_funcs_c.set('alg_res_func', c_body)

    def _make_dres_dc_function(self):
        py_body = []
        py_body.append('def dres_dc_function(time, dynamicVars, yprime, '
                       'constants):')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yprime = scipy.asarray(yprime)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('pd = scipy.zeros((%i, %i), scipy.float_)'
                       %(len(self.dynamicVars), len(self.dynamicVars)))
        py_body.append('')
        self._add_assignments_to_function_body(py_body)
        py_body.append('')

        c_body = []
        c_args = 'double *time_ptr, double *dynamicVars, double *yprime, '\
                'double *constants, double *pd'
        c_body.append('void dres_dc_function_(%s){' % c_args)
        self._prototypes_c.set('dres_dc_function', 
                               'void dres_dc_function_(%s);' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True)
        c_body.append('')

        N_dyn = len(self.dynamicVars)

        for res_ii, (rhs, dt) in enumerate(self._residual_terms):
            if rhs is None:
                rhs = '0'
            rhs_vars_used = ExprManip.extract_vars(rhs)
            for wrt_ii, wrt in enumerate(self.dynamicVars.keys()):
                pd_terms = []
                deriv = self.takeDerivative(rhs, wrt, rhs_vars_used)
                if (deriv != '0' and deriv != '0.0'):
                    py_body.append('# Residual %i wrt %s' % (res_ii, wrt))
                    py_body.append('pd[%i, %i] = %s' % (res_ii, wrt_ii, deriv))
                    deriv_c = ExprManip.make_c_compatible(deriv)
                    c_body.append('pd[%i] = %s;'
                                  % (res_ii+wrt_ii*N_dyn, deriv_c))

        py_body.append('')
        py_body.append('return pd')
        py_body = '\n    '.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('dres_dc_function', py_body)
        self._dynamic_funcs_c.set('dres_dc_function', c_body)

    def _make_dres_dcdot_function(self):
        py_body = []
        py_body.append('def dres_dcdot_function(time, dynamicVars, yprime, '
                       'constants):')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yprime = scipy.asarray(yprime)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('pd = scipy.zeros((%i, %i), scipy.float_)'
                       %(len(self.dynamicVars), len(self.dynamicVars)))
        py_body.append('')
        self._add_assignments_to_function_body(py_body)
        py_body.append('')

        c_body = []
        c_args = 'double *time_ptr, double *dynamicVars, double *yprime, '\
                'double *constants, double *pd'
        c_body.append('void dres_dcdot_function_(%s){' % c_args)
        self._prototypes_c.set('dres_dcdot_function', 
                               'void dres_dcdot_function_(%s);' % c_args)
        c_body.append('double time = *time_ptr;')
        c_body.append('')
        self._add_assignments_to_function_body(c_body, in_c=True)
        c_body.append('')

        N_dyn = len(self.dynamicVars)
        for res_ii, (rhs, dt) in enumerate(self._residual_terms):
            for wrt_ii, wrt in enumerate(self.dynamicVars.keys()):
                if dt == wrt:
                    py_body.append('# Derivative of residual term %i wrt '
                                   '%s_dot' % (res_ii, wrt))
                    py_body.append('pd[%i, %i] = -1' 
                                   % (res_ii, wrt_ii))
                    c_body.append('pd[%i] = -1;' % (res_ii+wrt_ii*N_dyn))

        py_body.append('')
        py_body.append('return pd')
        py_body = '\n    '.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('dres_dcdot_function', py_body)
        self._dynamic_funcs_c.set('dres_dcdot_function', c_body)

    def _make_ddaskr_jac(self):
        py_body = []
        py_body.append('def ddaskr_jac(time, dynamicVars, yprime, cj, '
                       'constants):')
        py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
        py_body.append('yprime = scipy.asarray(yprime)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('local_dres_dc = dres_dc_function(time, dynamicVars, '
                       'yprime, constants)')
        py_body.append('local_dres_dcdot = dres_dcdot_function(time, '
                       'dynamicVars, yprime, constants)')
        py_body.append('return local_dres_dc + cj * local_dres_dcdot')
        py_body = '\n    '.join(py_body)

        N_dyn = len(self.dynamicVars)
        c_body = []
        c_args = 'double *time_ptr, double *dynamicVars, double *yprime, '\
                'double *pd, double *cj_ptr, double *constants, int *intpar'
        c_body.append('void ddaskr_jac_(%s){' % c_args)
        self._prototypes_c.set('ddaskr_jac', 
                               'void ddaskr_jac_(%s);' % c_args)
        c_body.append('double cj = *cj_ptr;')
        c_body.append('')
        c_body.append('dres_dc_function_(time_ptr, dynamicVars, yprime, '
                      'constants, pd);')
        c_body.append('')
        # Like magic, this will initialize *all* elements to zero
        c_body.append('double local_dres_dcdot[%i*%i] = {0};' % (N_dyn, N_dyn))
        c_body.append('dres_dcdot_function_(time_ptr, dynamicVars, yprime, '
                      'constants, local_dres_dcdot);')
        c_body.append('')
        c_body.append('int ii;')
        c_body.append('for(ii=0; ii < %i; ii++){' % N_dyn**2)
        c_body.append('  pd[ii] += cj*local_dres_dcdot[ii];}')
        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('ddaskr_jac', py_body)
        self._dynamic_funcs_c.set('ddaskr_jac', c_body)

    def _make_dres_dsinglep(self):
        # loop over the optimizable params, and for each one create
        # its own dres_dparam function
        for wrt_ii, wrt in enumerate(self.optimizableVars.keys()):

            func_name = 'dres_d' + wrt

            py_body = []
            py_body.append('def ' + func_name + '(time, dynamicVars, '
                           'yprime, constants):')
            py_body.append('dynamicVars = scipy.asarray(dynamicVars)')
            py_body.append('constants = scipy.asarray(constants)')
            py_body.append('')
            py_body.append('pd = scipy.zeros(%i, scipy.float_)'
                           % len(self.dynamicVars))
            py_body.append('')

            self._add_assignments_to_function_body(py_body)
            py_body.append('')

            c_body = []
            c_args = 'double *time_ptr, double *dynamicVars, double *yprime, '\
                    'double *constants, double *pd'
            c_body.append('void ' + func_name + '_(%s){' % c_args)
            self._prototypes_c.set(func_name,
                                   'void ' + func_name + '_(%s);' % c_args)

            c_body.append('double time = *time_ptr;')
            c_body.append('')
            self._add_assignments_to_function_body(c_body, in_c=True)
            c_body.append('')

            N_dyn = len(self.dynamicVars)
            N_opt = len(self.optimizableVars)

            # Now the sensitivities. 
            # We'll cache lists of the variables in each rhs,
            # since extract_vars is slow.
            rhs_vars = {}
            for (rhs, dt) in self._residual_terms:
                if rhs is None:
                    rhs = '0'
                rhs_vars[rhs] = ExprManip.extract_vars(rhs)

            for res_ii, (rhs, dt) in enumerate(self._residual_terms):
                if rhs is None:
                    rhs = '0'
                rhs_vars_used = rhs_vars[rhs]
                deriv = self.takeDerivative(rhs, wrt, rhs_vars_used)
                if (deriv != '0' and deriv != '0.0'):
                    py_body.append('pd[%i] = %s' % (res_ii, deriv))
                    c_deriv = ExprManip.make_c_compatible(deriv)
                    c_body.append('pd[%i] = %s;' % (res_ii, c_deriv))
                              
            py_body.append('')
            py_body.append('return pd')
            py_body = '\n    '.join(py_body)

            c_body.append('}')
            c_body = os.linesep.join(c_body)

            self._dynamic_funcs_python.set(func_name, py_body)
            self._dynamic_funcs_c.set(func_name, c_body)

    def _make_sens_rhs(self):

        py_body = []
        py_body.append('def sens_rhs(time, sens_y, sens_yp, constants):')
        py_body.append('sens_y = scipy.asarray(sens_y)')
        py_body.append('sens_yp = scipy.asarray(sens_yp)')
        py_body.append('constants = scipy.asarray(constants)')
        py_body.append('')
        py_body.append('sens_res = scipy.zeros(%i, scipy.float_)'
                       % (2*len(self.dynamicVars)))
        py_body.append('')

        c_body = []
        c_args = 'double *time_ptr, double *sens_y, double *sens_yp, '\
                'double *cj_ptr, double *sens_res, int *ires_ptr, '\
                'double *constants, int *ipar'
        c_body.append('void sens_rhs_(%s){' % c_args)
        self._prototypes_c.set('sens_rhs', 
                               'void sens_rhs_(%s);' % c_args)
        c_body.append('')

        N_dyn = len(self.dynamicVars)
        N_const = len(self.constantVars)
        py_body.append('sens_res[:%(N_dyn)i] = res_function(time, '
                       'sens_y, sens_yp, constants)'
                       % {'N_dyn': N_dyn})
        py_body.append('')

        # for the python sens_rhs function, we add a dictionary
        # of dres_dparam function names indexed by parameter index.
        py_body.append('func_dict = {}')
        for wrt_ii, wrt in enumerate(self.optimizableVars.keys()):
            py_body.append('func_dict[' + str(wrt_ii) + '] = dres_d' + wrt)

        py_body.append('')

        py_body.append('p_index = int(constants[%i])' % N_const)

        # sens_rhs will access a particular dres_dparam function
        # based on the passed in parameter index
        py_body.append('dres_dp_func = func_dict[p_index]')

        py_body.append('')
        
        py_body.append('dc_dp = sens_y[%i:]' % N_dyn)
        py_body.append('dcdot_dp = sens_yp[%i:]' % N_dyn)
        py_body.append('local_dres_dp = dres_dp_func(time, sens_y, sens_yp, '
                       'constants)' % {'N_dyn': N_dyn, 'N_const': N_const})

        py_body.append('local_dres_dc = dres_dc_function(time, '
                       'sens_y, sens_yp, constants)' 
                       % {'N_dyn': N_dyn, 'N_const': N_const})
        py_body.append('local_dres_dcdot = dres_dcdot_function(time, '
                       'sens_y, sens_yp, constants)' 
                       % {'N_dyn': N_dyn, 'N_const': N_const})
        py_body.append('sens_res[%i:] = local_dres_dp '
                       '+ scipy.dot(local_dres_dc, dc_dp) '
                       '+ scipy.dot(local_dres_dcdot, dcdot_dp)' % N_dyn)

        # This fills in the first half of our sens_res
        c_body.append('res_function_(time_ptr, sens_y, sens_yp, cj_ptr, '
                      'sens_res, ires_ptr, constants, ipar);')
        c_body.append('')
        c_body.append('int p_index = (int)constants[%i];' % N_const)
        # we add an array that tracks the constants only, so that
        # the c version of the dres_dparam function gets an array of the
        # expected size
        c_body.append('double constants_only[%i];' % N_const)
        c_body.append('int jj;')
        c_body.append('for (jj = 0; jj < %s; jj++){' % N_const)
        c_body.append('constants_only[jj] = constants[jj];}')
        c_body.append('double *dc_dp = &sens_y[%i];' % N_dyn)
        c_body.append('double *dcdot_dp = &sens_yp[%i];' % N_dyn)
        # We'll directly fill dres_dp into the appropriate place in sens_res
        c_body.append('double *local_dres_dp = &sens_res[%i];' % N_dyn)
        c_body.append('int ii;')
        # sens_res isn't necessarily all zeros when passed in using the
        # cpointer.
        c_body.append('for(ii = 0; ii < %s; ii++){' % N_dyn)
        c_body.append('local_dres_dp[ii] = 0;}')

        # we add a switch-case statement to choose the proper dres_dparam
        # function based on the passed in parameted index
        c_body.append('switch(p_index)')
        c_body.append('{')
        for wrt_ii, wrt in enumerate(self.optimizableVars.keys()):
            c_body.append('case ' + str(wrt_ii) + ' : dres_d' + wrt + '_(time_ptr, sens_y, sens_yp, constants_only, local_dres_dp);')
            c_body.append('break;'), 
        c_body.append('}')
        # Fill in dres_dc
        c_body.append('double local_dres_dc[%i] = {0};' % N_dyn**2)
        c_body.append('dres_dc_function_(time_ptr, sens_y, sens_yp, constants, '
                      'local_dres_dc);')
        c_body.append('int row, col;')
        c_body.append('for(row = 0; row < %i; row++){' % N_dyn)
        c_body.append('for(col = 0; col < %i; col++){' % N_dyn)
        c_body.append('sens_res[row+%i] += local_dres_dc[row + col*%i]*dc_dp[col];}}' % (N_dyn, N_dyn))

        c_body.append('double local_dres_dcdot[%i] = {0};' % N_dyn**2)
        c_body.append('dres_dcdot_function_(time_ptr, sens_y, sens_yp, '
                      'constants, local_dres_dcdot);')
        c_body.append('for(row = 0; row < %i; row++){' % N_dyn)
        c_body.append('for(col = 0; col < %i; col++){' % N_dyn)
        c_body.append('sens_res[row+%i] += local_dres_dcdot[row + col*%i]'
                      '*dcdot_dp[col];}}' % (N_dyn, N_dyn))

        py_body.append('')
        py_body.append('return sens_res')
        py_body = '\n    '.join(py_body)

        c_body.append('}')
        c_body = os.linesep.join(c_body)

        self._dynamic_funcs_python.set('sens_rhs', py_body)
        self._dynamic_funcs_c.set('sens_rhs', c_body)

    def _make_log_funcs(self):
        N_dyn = len(self.dynamicVars)

        py_body = []
        py_body.append('def res_function_logdv(time, log_dv, log_yp, '
                       'constants):')
        py_body.append('log_dv = scipy.asarray(log_dv)')
        py_body.append('log_yp = scipy.asarray(log_yp)')
        py_body.append('dynamicVars = scipy.exp(log_dv)')
        py_body.append('dynamicVars = scipy.maximum(dynamicVars, %g)' 
                       % _double_tiny_)
        py_body.append('yprime = log_yp * dynamicVars')
        py_body.append('return res_function(time, dynamicVars, yprime, '
                       'constants)')
        py_body = '\n    '.join(py_body)
        self._dynamic_funcs_python.set('res_function_logdv', py_body)

        c_body = []
        c_args = 'double *time_ptr, double *log_dv, double *log_yp, '\
                'double *cj_ptr, double *residual, int *ires_ptr, '\
                'double *constants, int *ipar'
        c_body.append('void res_function_logdv_(%s){' % c_args)
        self._prototypes_c.set('res_function_logdv',
                               'void res_function_logdv_(%s);' % c_args)
        c_body.append('double dynamicVars[%i];' % N_dyn)
        c_body.append('double yprime[%i];' % N_dyn)
        c_body.append('int ii;')
        c_body.append('for(ii = 0; ii < %i; ii++){' % N_dyn)
        c_body.append('dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);')
        c_body.append('yprime[ii] = log_yp[ii] * dynamicVars[ii];}')
        c_body.append('res_function_(time_ptr, dynamicVars, yprime, cj_ptr, '
                      'residual, ires_ptr, constants, ipar);')
        c_body.append('}')
        c_body = os.linesep.join(c_body)
        self._dynamic_funcs_c.set('res_function_logdv', c_body)

        py_body = []
        py_body.append('def root_func_logdv(time, log_dv, log_yp, constants):')
        py_body.append('log_dv = scipy.asarray(log_dv)')
        py_body.append('log_yp = scipy.asarray(log_yp)')
        py_body.append('dynamicVars = scipy.exp(log_dv)')
        py_body.append('dynamicVars = scipy.maximum(dynamicVars, %g)'
                       % _double_tiny_)
        py_body.append('yprime = log_yp * dynamicVars')
        py_body.append('return root_func(time, dynamicVars, yprime, '
                       'constants)')
        py_body = '\n    '.join(py_body)
        self._dynamic_funcs_python.set('root_func_logdv', py_body)

        c_body = []
        c_args = 'int *neq_ptr, double *time_ptr, double *log_dv, '\
                'double *log_yp, int *nrt_ptr, double *root_devs, '\
                'double *constants, int *ipar'
        c_body.append('void root_func_logdv_(%s){' % c_args)
        self._prototypes_c.set('root_func_logdv',
                               'void root_func_logdv_(%s);' % c_args)
        c_body.append('double dynamicVars[%i];' % N_dyn)
        c_body.append('double yprime[%i];' % N_dyn)
        c_body.append('int ii;')
        c_body.append('for(ii = 0; ii < %i; ii++){' % N_dyn)
        c_body.append('dynamicVars[ii] = max(exp(log_dv[ii]), DBL_MIN);')
        c_body.append('yprime[ii] = log_yp[ii] * dynamicVars[ii];}')
        c_body.append('root_func_(neq_ptr, time_ptr, dynamicVars, yprime, '
                      'nrt_ptr, root_devs, constants, ipar);')
        c_body.append('}')
        c_body = os.linesep.join(c_body)
        self._dynamic_funcs_c.set('root_func_logdv', c_body)


        py_body = []
        py_body.append('def sens_rhs_logdv(time, sens_y, sens_yp, constants):')
        # We need to copy sens_y and sens_yp here, since we aren't allowed to
        # change their values.
        py_body.append('sens_y = scipy.array(sens_y)')
        py_body.append('sens_yp = scipy.array(sens_yp)')
        py_body.append('sens_y[:%i] = scipy.exp(sens_y[:%i])' % (N_dyn, N_dyn))
        py_body.append('sens_y[:%i] = scipy.maximum(sens_y[:%i], %g)'
                       % (N_dyn, N_dyn, _double_tiny_))
        py_body.append('sens_yp[:%i] = sens_yp[:%i] * sens_y[:%i]'
                       % (N_dyn, N_dyn, N_dyn))
        py_body.append('return sens_rhs(time, sens_y, sens_yp, constants)')
        py_body = '\n    '.join(py_body)
        self._dynamic_funcs_python.set('sens_rhs_logdv', py_body)

        c_body = []
        c_args = 'double *time_ptr, double *sens_y_log, double *sens_yp_log, '\
                'double *cj_ptr, double *sens_res, int *ires_ptr, '\
                'double *constants, int *ipar'
        c_body.append('void sens_rhs_logdv_(%s){' % c_args)
        self._prototypes_c.set('sens_rhs_logdv_',
                               'void sens_rhs_logdv_(%s);' % c_args)
        c_body.append('double sens_y[%i];' % (N_dyn*2))
        c_body.append('double sens_yp[%i];' % (N_dyn*2))
        c_body.append('int ii;')
        c_body.append('for(ii = 0; ii < %i; ii++){' % N_dyn)
        c_body.append('sens_y[ii] = max(exp(sens_y_log[ii]), DBL_MIN);')
        c_body.append('sens_yp[ii] = sens_yp_log[ii] * sens_y[ii];}')
        c_body.append('for(ii = %i; ii < %i; ii++){' % (N_dyn, 2*N_dyn))
        c_body.append('sens_y[ii] = sens_y_log[ii];')
        c_body.append('sens_yp[ii] = sens_yp_log[ii];}')
        c_body.append('sens_rhs_(time_ptr, sens_y, sens_yp, cj_ptr, sens_res, '
                      'ires_ptr, constants, ipar);')
        c_body.append('}')
        c_body = os.linesep.join(c_body)
        self._dynamic_funcs_c.set('sens_rhs_logdv', c_body)

    _bypass_integrate_stochastic_tidbit = False
    def _make_integrate_stochastic_tidbit(self):
        # Function definitions
        py_body = []
        py_body.append('def integrate_stochastic_tidbit'
                       '(seed, reseed, time, dv, cv, rmsd, stop_time):')

        c_args = 'unsigned long* seed_ptr, int* reseed, double* time_ptr, '\
                 'int* dv, double* cv, '\
                 'double* rmsd_ptr, double* stop_time_ptr, double* trajectory'
        c_body = []
        c_body.append('void integrate_stochastic_tidbit_(%s) {'%c_args)

        if self._bypass_integrate_stochastic_tidbit:
            py_body.append('pass')
            py_body = '\n    '.join(py_body)

            c_body.append('}')
            c_body = os.linesep.join(c_body)

            self._prototypes_c.set('integrate_stochastic_tidbit', 
                                   'void integrate_stochatic_tidbit_(%s);' 
                                   % c_args)

            self._dynamic_funcs_python.set('integrate_stochastic_tidbit', 
                                           py_body)
            self._dynamic_funcs_c.set('integrate_stochastic_tidbit', c_body)
            return 

        # Test for conditions under which we cannot make these functions
        err_body = []
        if len(self.reactions) == 0:
            err_body.append('# No reactions present.')
        if len(self.rateRules) > 0:
            err_body.append('# %i rate rules are present'%len(self.rateRules))
        if len([_ for _ in self.events if not _.timeTriggered]) > 0:
            err_body.append('# %i non-time-triggered events are present'%\
                            len([_ for _ in self.events \
                                 if not _.timeTriggered]))
        if len(self.algebraicRules) > 0:
            err_body.append('# %i algebraic rules are present'%
                            len(self.algebraicRules))
        if len(self.algebraicVars) > 0:
            err_body.append('# %i algebraic variables are present'%
                            len(self.algebraicRules))
        for rxn in self.reactions.values():
            for val in rxn.stoichiometry.values():
                try: intVal = int(val)
                except: intVal = None
                if intVal is None:
                    err_body.append('# Noncastable stoichiometry for %s: %s'%
                                    (rxn.id, repr(val)))
                elif intVal != val:
                    err_body.append('# Reaction stoichiometry for %s: %s!=%s'%
                                    (rxn.id, repr(val), repr(intVal)))

        if err_body!=[]:
            py_body = py_body + err_body
            py_body.append('pass')
            py_body = '\n    '.join(py_body)
            self._dynamic_funcs_python.set('integrate_stochastic_tidbit',
                                           py_body)
            
            c_body.append('return;}')
            c_body = os.linesep.join(c_body)
            self._prototypes_c.set('integrate_stochastic_tidbit', 
                                   'void integrate_stochatic_tidbit_(%s);' %
                                   c_args)
            self._dynamic_funcs_c.set('integrate_stochastic_tidbit', c_body)

            return
        
        N_CON = len(self.constantVars)
        N_DYN = len(self.dynamicVars)
        N_ASN = len(self.assignedVars)
        N_RXN = len(self.reactions)

        # Use a custom AST walk to get 12x faster than ExprManip.sub_for_vars
        mapping = {}
        for index, name in enumerate(self.constantVars.keys()):
            mapping[name] = 'cv[%i]'%index
        for index, name in enumerate(self.dynamicVars.keys()):
            mapping[name] = 'dv[%i]'%index
        for index, name in enumerate(self.assignedVars.keys()):
            mapping[name] = 'av[%i]'%index
        
        class Parse(ExprManip.AST.compiler.visitor.ASTVisitor):
            def __call__(slf,s,c=True):
                if c: s = ExprManip.make_c_compatible(s)
                ast = ExprManip.strip_parse(s)
                ExprManip.AST.compiler.walk(ast, slf)
                return ExprManip.ast2str(ast)

            def visitName(slf,node,*args):
                if mapping.has_key(node.name): node.name=mapping[node.name]
        parse = Parse()
        
        # Find the stoichiometry and reaction dependencies
        asnKeys = self.assignedVars.keys()
        evars = [] # The extracted dynamic variables each rxn depends upon
        for rxnInd, rxn in enumerate(self.reactions.values()):
            evar = ExprManip.Extraction.extract_vars(rxn.kineticLaw)
            evar = sets.Set(evar)
            while len(evar.intersection(asnKeys)):
                for asnName in evar.intersection(asnKeys):
                    rule = self.assignmentRules.get(asnName)
                    evar.remove(asnName)
                    evar.union_update( \
                        ExprManip.Extraction.extract_vars(rule))
            evars.append(evar)

        stch_body = []
        depd_body = []
        for rxnInd, (rxn_id, rxn) in enumerate(self.reactions.items()):
            stch_body.append(repr([int(rxn.stoichiometry.get(_,0))
                                   for _ in self.dynamicVars.keys()])[1:-1])
            varNames = sets.Set([n
                                 for n,v in rxn.stoichiometry.items() if v!=0])
            depd = [int(len(evar.intersection(varNames))>0) for evar in evars]
            depd_body.append(repr(depd)[1:-1])
        if N_RXN>0: depd_body.append(repr([1]*N_RXN)[1:-1])

        # Python body assembly
        py_body.append('import numpy.random')
        py_body.append('if reseed: rs = numpy.random.seed(seed)')
        py_body.append('')
        if N_RXN > 0:
            py_body.append('stch = [[%s]]' % '],\n        ['.join(stch_body))
            py_body.append('depd = [[%s]]' % '],\n        ['.join(depd_body))
        else:
            py_body.append('stch = [[]]')
            py_body.append('depd = [[]]')
        py_body.append('')
        py_body.append('av = [0.0]*%i'%N_ASN)
        py_body.append('dv0 = [_ for _ in dv]')
        py_body.append('props = [0.0]*%i'%N_RXN)
        py_body.append('rxnInd = %i'%N_RXN)
        py_body.append('while time < stop_time:')
        for index, rule in enumerate(self.assignmentRules.values()):
            py_body.append('    av[%i]=%s;'%(index, parse(rule, False)))
        py_body.append('')
        for index, rxn in enumerate(self.reactions.values()):
            py_body.append('    if depd[rxnInd][%i]: props[%i]=%s'%
                           (index, index, parse(rxn.kineticLaw, False)))
        py_body.append('    propensity=sum(props)')
        py_body.append('    if propensity==0.0:')
        py_body.append('        dt = stop_time - time')
        py_body.append('        time=stop_time')
        py_body.append('        break')
        py_body.append('    ')
        py_body.append('    dt = -log(1.0-numpy.random.rand())/propensity')
        py_body.append('    time += dt')
        py_body.append('    ')
        py_body.append('    selection=propensity*numpy.random.rand()')
        py_body.append('    for rxnInd in range(%i):'%N_RXN)
        py_body.append('        if selection<props[rxnInd]: break')
        py_body.append('        else: selection-=props[rxnInd]')
        py_body.append('    ')
        py_body.append('    for index, st in enumerate(stch[rxnInd]):')
        py_body.append('        dv[index] += st')
        py_body.append('    ')
        py_body.append('    _rmsd = sum([(_-__)**2 for _, __ in zip(dv,dv0)])')
        py_body.append('    _rmsd = (_rmsd/%i.)**0.5'%N_DYN)
        py_body.append('    if _rmsd>rmsd: break')
        py_body.append('')
        py_body.append('if time > stop_time: # Back interpolate')
        py_body.append('    trajectory = ['
                       'float(_)+float(stch[rxnInd][ii])/dt*(stop_time-time) '
                       'for ii, _ in enumerate(dv)]')
        py_body.append('else:')
        py_body.append('    stop_time = time')
        py_body.append('    trajectory = [float(_) for _ in dv]')
        py_body.append('return time, dv, stop_time, trajectory')
        py_body = '\n    '.join(py_body)

        # C body assembly
        c_body.append('  int i; /* Temp variables */')
        c_body.append('')
        c_body.append('  unsigned long seed = *seed_ptr;')
        c_body.append('  if (*reseed) {init_genrand(seed);}')
        c_body.append('')
        if N_RXN > 0:
            c_body.append('  short stch[%i][%i] = {{%s}};' %
                          (N_RXN, N_DYN,
                           '},\n                    {'.join(stch_body)))
            c_body.append('  short depd[%i+1][%i] = {{%s}};' %
                          (N_RXN, N_RXN,
                           '},\n                    {'.join(depd_body)))
        else:
            c_body.append('  short stch[0][0];')
        c_body.append('')
        c_body.append('  double time = *time_ptr;')
        c_body.append('  double sd = (*rmsd_ptr)*(*rmsd_ptr)*%i.;'%N_DYN)
        c_body.append('  double stop_time = *stop_time_ptr;')
        c_body.append('  double dt=0.0;')
        c_body.append('')
        c_body.append('  double dv0[%i];'%N_DYN)
        c_body.append('  for (i=0;i<%i;i++) {dv0[i]=dv[i];}'%N_DYN)
        c_body.append('')
        c_body.append('  int rxnInd = %i;'%N_RXN)
        c_body.append('  double propensity, selection, props[%i], av[%i];'%
                      (N_RXN, N_ASN))
        c_body.append('  while (time < stop_time) {')
        for index, rule in enumerate(self.assignmentRules.values()):
            c_body.append('    av[%i]=%s;'%(index, parse(rule)))
        c_body.append('')
        for index, rxn in enumerate(self.reactions.values()):
            c_body.append("    if (depd[rxnInd][%i]) {props[%i]=%s;}"% \
                          (index, index, parse(rxn.kineticLaw)))
        c_body.append('')
        c_body.append('    propensity = 0.0;')
        c_body.append('    for (i=0;i<%i;i++) {'%N_RXN)
        c_body.append('      propensity += props[i];}')
        c_body.append('   if (propensity<=0.0) {')
        c_body.append('      dt = stop_time-time;')
        c_body.append('      time = stop_time;')
        c_body.append('      break;')
        c_body.append('   }')
        c_body.append('')
        c_body.append('    dt = -log(1.0-genrand_real32())/propensity;')
        c_body.append('    time += dt;')
        c_body.append('')
        c_body.append('    selection = propensity * genrand_real32();')
        c_body.append('')
        c_body.append('    for (rxnInd=0; rxnInd<%i; rxnInd++) {'%N_RXN)
        c_body.append('      if (selection < props[rxnInd]) {break;}')
        c_body.append('      else {selection -= props[rxnInd];}}')
        c_body.append('')
        c_body.append('    for (i=0;i<%i;i++) {dv[i]+=stch[rxnInd][i];}'%N_DYN)
        c_body.append('')
        c_body.append('    double _sd = 0.0;')
        c_body.append('    for (i=0;i<%i;i++) {'%N_DYN)
        c_body.append('        _sd += (dv0[i]-dv[i])*(dv0[i]-dv[i]);')
        c_body.append('    }')
        c_body.append('    if (_sd > sd) {break;}')
        c_body.append('  }')
        c_body.append('')
        c_body.append('  for (i=0;i<%i;i++) {'%N_DYN)
        c_body.append('    trajectory[i]=(double)dv[i];')
        c_body.append('    if (time > stop_time) {')
        c_body.append('      trajectory[i] += (double)stch[rxnInd][i]/dt*'
                      '(stop_time - time);')
        c_body.append('    }')
        c_body.append('  }')
        c_body.append('  if (time>stop_time) {(*stop_time_ptr) = stop_time;}')
        c_body.append('  else {(*stop_time_ptr) = time;}')
        c_body.append('')
        c_body.append('  (*time_ptr) = time;')
        c_body.append('}')
        c_body.append('')
        c_body = os.linesep.join(c_body)

        self._prototypes_c.set('integrate_stochastic_tidbit', 
                               'void integrate_stochatic_tidbit_(%s);' % c_args)

        self._dynamic_funcs_python.set('integrate_stochastic_tidbit', py_body)
        self._dynamic_funcs_c.set('integrate_stochastic_tidbit', c_body)


    def _add_assignments_to_function_body(self, body, in_c = False,
                                          include_dts=False):
        """
        Adds the assignment rules for this Network to the list of lines that
        are in body.
        """
        # We make sure these guys are arrays, so that we can be confident
        #  in using .item().
        # We loop to assign our constantVarValues and our dynamicVars
        #  to local variables for speed (avoid repeated accesses) and
        #  for readability
        for arg, var_names in zip(['constants', 'dynamicVars'],
                                  [self.constantVars.keys(),
                                   self.dynamicVars.keys()]):
            for ii, id in enumerate(var_names):
                if not in_c:
                    body.append('%s = %s.item(%i)' % (id, arg, ii))
                    if include_dts and arg != 'constants':
                        body.append('%s_deriv_wrt_time = yprime.item(%i)' 
                                    % (id,ii))
                else:
                    body.append('double %s = %s[%i];' % (id, arg, ii))
                    if include_dts and arg != 'constants':
                        body.append('double %s_deriv_wrt_time = yprime[%i];' 
                                    % (id,ii))
            body.append('')

        for variable, math in self.assignmentRules.items():
            if not in_c:
                body.append('%s = %s' % (variable, math))
            else:
                c_math = ExprManip.make_c_compatible(math)
                body.append('double %s = %s;' % (variable, c_math))
            if include_dts:
                # This zero is in case assigment doesn't depend on any
                # dynamic variables. In that case, simply_expr could fail
                # on an empty string.
                rhs_terms = ['0']
                # Which variables does assignment rule depend on?
                dep_vars = ExprManip.extract_vars(math)

                # Do the chain rule for dynamic variables.
                for dynvar in dep_vars.intersection(self.dynamicVars.keys()):
                    dmath_dvar = ExprManip.diff_expr(math, dynvar)
                    rhs_terms.append('%s * %s_deriv_wrt_time'
                                     % (dmath_dvar, dynvar))
                # Add in explicit time dependence.
                if 'time' in dep_vars:
                    dmath_dtime = ExprManip.diff_expr(math, 'time')
                    rhs_terms.append(dmath_dtime)

                rhs = ' + '.join(rhs_terms)
                rhs = ExprManip.simplify_expr(rhs)
                if not in_c:
                    body.append('%s_deriv_wrt_time = %s' % (variable, rhs))
                else:
                    c_rhs = ExprManip.make_c_compatible(rhs)
                    body.append('double %s_deriv_wrt_time = %s;'
                                % (variable, c_rhs))

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
        var_vals = [(id, self.get_var_val(id)) for id in self.variables.keys()]
        var_vals = dict(var_vals)
        var_vals['time'] = time_fired
        # Calculate the values that will get assigned to the variables
        for id, rhs in event.event_assignments.items():
            new_values[id] = self.evaluate_expr(rhs, time_fired, var_vals)

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

    # We cache derivatives wrt dynamicVariables, because they'll be redundantly
    # taken many, many times otherwise, and it can be a major hold-up.
    # Unforunately, the structures we want to index by aren't hashable, so this
    # can't be a dictionary. Luckily that shouldn't slow things down much,
    # because we shouldn't have to sort through many structures in our
    # applications.
    _sens_event_derivs_cache = []
    def executeEventAndUpdateSens(self, holder, ysens_pre_exec, opt_var):
        curr_struct = self._get_structure()
        curr_events = copy.deepcopy(self.events)
        # Search through the _sens_event_derivs_cache for the derivs
        # corresponding to this structure.
        for struct, events, derivs_cache in self._sens_event_derivs_cache:
            if curr_struct == struct and curr_events == events:
                break
        else:
            derivs_cache = {}
            self._sens_event_derivs_cache.append((curr_struct, curr_events,
                                                 derivs_cache))

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
        holder.ysens_pre_exec = ysens_pre_exec

        if not event.event_assignments:
            # Event does nothing, so it can't affect parameter senstivities
            # Note that we don't get derivatives of time of firing and execution
            #  if we skip this. However, it is very useful for when we have
            #  events monitoring derivatives of variables, because then we
            #  can't calculate those event time derivatives.
            holder.ysens_post_exec = ysens_pre_exec
            return holder.ysens_post_exec

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
        var_vals_fired = [(id, self.get_var_val(id)) for id 
                          in self.variables.keys()]
        var_vals_fired = dict(var_vals_fired)
        var_vals_fired['time'] = time
        # We need to do this if we don't have a chained event.
        if dtf_dp is None:
            # We convert our trigger to smooth function we can take the 
            # derivative of so we can calculate sensitivity of firing time.
            # The trigger is passed to us because it might be a clause of a
            # more complicated trigger.
            trigger = self.event_clauses[clause_index]
            trigger = self._smooth_trigger(trigger)
            trigger_vars_used = ExprManip.extract_vars(trigger)

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
                try:
                    dtrigger_ddv = derivs_cache[(trigger, id)]
                except KeyError:
                    dtrigger_ddv = self.takeDerivative(trigger, id, 
                                                       trigger_vars_used)
                    derivs_cache[(trigger, id)] = dtrigger_ddv
                dtrigger_ddv = self.evaluate_expr(dtrigger_ddv, time_fired,
                                                  var_vals_fired)
                # This is the index of ddv_dp in the y array
                index_in_y = ii + N_dv
                ddv_dp = ysens_fired[index_in_y]
                ddv_dt = yp_fired[ii]
                numerator += dtrigger_ddv * ddv_dp
                denominator += dtrigger_ddv * ddv_dt

            # Second term in numerator: dtrigger_dp
            dtrigger_dp = self.takeDerivative(trigger, opt_var, 
                                              trigger_vars_used)
            dtrigger_dp = self.evaluate_expr(dtrigger_dp, time_fired, 
                                             var_vals_fired)
            numerator += dtrigger_dp

            # Second term in denominator: dtrigger_dt
            dtrigger_dt = self.takeDerivative(trigger, 'time', 
                                              trigger_vars_used)
            dtrigger_dt = self.evaluate_expr(dtrigger_dt, time_fired,
                                             var_vals_fired)
            denominator += dtrigger_dt

            dtf_dp = -numerator/denominator

        # Now compute the sensitivity of the delay
        #
        # ddelay/dp = ddelay/dp + ddelay_dy*dy/dp + ddelay_dy*dy/dt*dtf/dp
        #             + ddelay/dt*dtf/dp
        delay = event.delay
        if not isinstance(delay, str):
            ddelay_dp = 0
        else:
            # First term
            delay_vars_used = ExprManip.extract_vars(delay)
            ddelay_dp = self.takeDerivative(delay, opt_var, delay_vars_used)
            ddelay_dp = self.evaluate_expr(ddelay_dp, time_fired,
                                           var_vals_fired)

            # Second and third terms: ddelay_dy*dy/dp + ddelay_dy*dy/dt*dtf/dp
            for ii, id in enumerate(self.dynamicVars.keys()):
                try:
                    ddelay_dy = derivs_cache[(delay, id)]
                except KeyError:
                    ddelay_dy = self.takeDerivative(delay, id, delay_vars_used)
                    derivs_cache[(delay, id)] = ddelay_dy
                ddelay_dy = self.evaluate_expr(ddelay_dy, time_fired,
                                               var_vals_fired)
                # This is the index of ddv_dp in the y array
                index_in_y = ii + N_dv
                dy_dp = ysens_fired[index_in_y]
                dy_dt = yp_fired[ii]
                ddelay_dp += ddelay_dy * dy_dp
                ddelay_dp += ddelay_dy * dy_dt * dtf_dp

            # Fourth term: ddelay/dt*dtf/dp
            ddelay_dt = self.takeDerivative(delay, 'time', delay_vars_used)
            ddelay_dt = self.evaluate_expr(ddelay_dt, time_fired, 
                                           var_vals_fired)
            ddelay_dp += ddelay_dt * dtf_dp

        # This is the delay in the time of execution
        dte_dp = dtf_dp + ddelay_dp

        # We store d(time of execution) for use by other events that may be
        #  chained to this one.
        holder.dte_dp = dte_dp

        # We update our Network so we can get all the variable values
        #  that were just prior to the event executing.
        self.updateVariablesFromDynamicVars(ysens_pre_exec, time_exec)
        var_vals_pre_exec = [(id, self.get_var_val(id)) for id 
                          in self.variables.keys()]
        var_vals_pre_exec = dict(var_vals_pre_exec)
        var_vals_pre_exec['time'] = time
        
        # Now compute the sensitivity of each of our new values
        for y_ii, y_id in enumerate(self.dynamicVars.keys()):
            if not event.event_assignments.has_key(y_id):
                # This is the index of y's sensitivity in the sensitivity array
                index_in_y = y_ii + N_dv
                # dy_dp after the event of course begins equal to the current
                #  sensitivity
                dy_dp_post_exec = ysens_pre_exec[index_in_y]
                # These next two terms account for possible changes in the 
                #  derivative of y due to event execution.
                # dy/dt(pre_exec) * dte/dp term
                dy_dp_post_exec += yp_pre_exec[y_ii] * dte_dp
                # -dy/dt(post_exec) * dte/dp
                dy_dp_post_exec -= yp_post_exec[y_ii] * dte_dp
                ysens_post_exec[index_in_y] = dy_dp_post_exec
            else:
                # If y has been assigned the value a, then the sensitivity of 
                #  y after the event is da/dp - dy/dt(post_exec) * dte/dp
                a = event.event_assignments.get(y_id)
                if not isinstance(a, str):
                    # If a isn't a string, it must be a constant, so da/dp = 0.
                    dy_dp_post_exec = -yp_post_exec[y_ii] * dte_dp
                    index_in_y = y_ii + N_dv
                    ysens_post_exec[index_in_y] = dy_dp_post_exec
                else:
                    # da/dp= da/dp + da/dy * dy/dp + da/dy * dy/dt * dtf/dp
                    #                   + da/dt * dtf/dp
                    # All these are calculated at the firing time
            
                    # We only need to deal with whatever clause in the 
                    # piecewise is applicable.
                    a = self._sub_for_piecewise(a, time_fired)
                    a_vars_used = ExprManip.extract_vars(a)

                    # First term: da/dp 
                    da_dp = self.takeDerivative(a, opt_var, a_vars_used)
                    da_dp = self.evaluate_expr(da_dp, time_fired, 
                                               var_vals_fired)
                    dy_dp_post_exec = da_dp

                    # Second term: da/dy * dy/dp 
                    # Third term: da/dy * dy/dt * dtf/dp 
                    # Together = da/dy * (dy/dp + dy/dt * dtf/dp)
                    for other_ii,other_id in enumerate(self.dynamicVars.keys()):
                        try:
                            da_dother = derivs_cache[(a, other_id)]
                        except KeyError:
                            da_dother = self.takeDerivative(a, other_id, 
                                                            a_vars_used)
                            derivs_cache[(a, other_id)] = da_dother
                        da_dother = self.evaluate_expr(da_dother, time_fired, 
                                                       var_vals_fired)

                        # This is the index of dother_dp in the y array
                        index_in_y = other_ii + N_dv
                        dother_dp = ysens_fired[index_in_y]

                        dother_dt = yp_fired[other_ii]

                        dy_dp_post_exec += da_dother * dother_dp
                        dy_dp_post_exec += da_dother * dother_dt * dtf_dp

                    # Fourth term: da/dt * dtf/dp
                    da_dt = self.takeDerivative(a, 'time', a_vars_used)
                    da_dt = self.evaluate_expr(da_dt, time_fired, 
                                               var_vals_fired)
                    dy_dp_post_exec += da_dt * dtf_dp

                    # Finally, -dy/dt(post_exec) * dte/dp
                    # We calculated dy/dt(new) up above in this function
                    dy_dp_post_exec -= yp_post_exec[y_ii] * dte_dp

                    index_in_y = y_ii + N_dv
                    ysens_post_exec[index_in_y] = dy_dp_post_exec

        holder.ysens_post_exec = ysens_post_exec
        return ysens_post_exec

    def _manualCrossReferences(self, flag=True):
        """
        Disable or enable automatic evaluation of cross-reference lists.

        Delaying cross reference creation can dramatically speed up creation
        of complex Networks. Once all Network components have been added,
        net._makeCrossReferences must be called manually. You may also need
        to call updateAssignedVars.

        _makeCrossReferences must also be called if the constancy or
        optimizability of a variable is changed.

        This is an advanced option, be careful!
        """
        self._manualCrossReferences_flag=flag

    #
    # Internally useful things
    #
    def _makeCrossReferences(self):
        """
        Create the cross-reference lists for the Network.
        """

        # Once makeCrossReferences is explicitly called, we assume that it
        # should be automatically called
        self._manualCrossReferences(flag=False)
        
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

        # Set the algebraicVars list to the list we just compiled, after sorting
        #  based on the order in self.dynamicVars
        sorted_alg_vars = [(id, self.variables.get(id))
                           for id in self.dynamicVars.keys()
                           if id in vars_in_alg_rules]
        self.algebraicVars = KeyedList(sorted_alg_vars)


    def disable_deriv_funcs(self):
        """
        Disables the creation of derivative functions, which can speed up
        integration.
        """
        
        if self.deriv_funcs_enabled == True:
            self.ddaskr_jac = None
            self._dynamic_structure_methods = ['_make_res_function',
                                          '_make_alg_deriv_func',
                                          '_make_alg_res_func',
                                          '_make_integrate_stochastic_tidbit'
                                          ]
            self.deriv_funcs_enabled = False

    def enable_deriv_funcs(self):
        """
        Enables the creation of derivative functions, which allows
        sensitivity integration and may speed up normal integration
        as well.
        """        

        if self.deriv_funcs_enabled == False:
            self._dynamic_structure_methods = ['_make_res_function',
                                          '_make_alg_deriv_func',
                                          '_make_alg_res_func',
                                          '_make_dres_dc_function',
                                          '_make_dres_dcdot_function',
                                          '_make_ddaskr_jac',
                                          '_make_dres_dsinglep',
                                          '_make_sens_rhs',
                                          '_make_log_funcs',
                                          '_make_integrate_stochastic_tidbit'
                                          ]
            self.deriv_funcs_enabled = True
    

    _last_structure = None
    # This is an option to disable compilation of C modules.
    disable_c = SloppyCell.disable_c
    _last_disabled_c = disable_c
    def compile(self, disable_c=None):
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
        if disable_c is None:
            disable_c = self.disable_c

        curr_structure = self._get_structure()
        structure_changed = (curr_structure != self._last_structure)
        if structure_changed:
            self._last_structure = copy.deepcopy(curr_structure)
            # The dynamic function lists need to be reset because the number
            #  and names of the dres_dparam functions can differ between
            #  networks.
            # We do not include the list of function definitions because
            #  they are handled differently than other dynamic functions.
            # Users should make sure that if they remove a function definition
            #  from a compiled network that they also remove references to
            #  that function from the rest of the network

            # We clear out all the dynamic functions that have been defined.
            all_dynamic_keys = sets.Set(self._dynamic_funcs_python.keys())
            all_dynamic_keys.union_update(self._dynamic_funcs_c.keys())
            for dynamic_func in all_dynamic_keys:
                try:
                    delattr(self, dynamic_func)
                except AttributeError:
                    pass

            self._dynamic_funcs_python = KeyedList()
            self._prototypes_c = self._common_prototypes_c.copy()
            self._dynamic_funcs_c = KeyedList()

        reexec = False
        if structure_changed:
            logger.debug('Network %s: compiling structure.' % self.id)
            self._makeCrossReferences()
            # Check whether system appears well-determined.
            if len(self.algebraicVars) > len(self.algebraicRules):
                raise ValueError('System appears under-determined. '
                                 'Not enough  algebraic rules. '
                                 '\n algebraicVars (%i) = \n%s '
                                 '\n algebraicRules (%i) = \n%s '%(len(self.algebraicVars), self.algebraicVars,
                                                                   len(self.algebraicRules), self.algebraicRules)) 
            self.make_func_defs()
            self._makeDiffEqRHS()
            for method in self._dynamic_structure_methods:
                getattr(self, method)()
            reexec = True

        if self.events != getattr(self, '_last_events', None)\
           or structure_changed:
            logger.debug('Network %s: compiling events.' % self.id)
            self._makeCrossReferences()
            for method in self._dynamic_event_methods: 
                getattr(self, method)()
            self._last_events = copy.deepcopy(self.events)
            reexec = True

        # after compile, check that there are not more algebraic variables than
        # algebraic rules.
        if len(self.algebraicVars) > len(self.algebraicRules):
            raise ValueError('System appears under-determined. '
                             'Not enough  algebraic rules. '
                             '\n algebraicVars (%i) = \n%s '
                             '\n algebraicRules (%i) = \n%s '%(len(self.algebraicVars), self.algebraicVars,
                                                               len(self.algebraicRules), self.algebraicRules)) 

        if self._last_disabled_c != disable_c:
            self._last_disabled_c = disable_c
            reexec = True

        self.compiled=True
        if reexec:
            self.exec_dynamic_functions(disable_c=disable_c)

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
                     self.algebraicRules, self.deriv_funcs_enabled)
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

    # We cache our exec'd python functions and compiled C modules to minimize
    #  redundant compiles.
    _py_func_dict_cache = {}
    _c_module_cache = {}
    # This is an option to disable compilation of C modules.
    def exec_dynamic_functions(self, disable_c=False, del_c_files=True, 
                               curr_c_code=None):
        # only get the bodies that were created.
        curr_py_bodies = '\n'.join([body for body in self._dynamic_funcs_python.values()\
                                    if body != None])

        key = (curr_py_bodies, tuple(self._func_strs.items()))
        # Search our cache of python functions.
        try:
            py_func_dict = self._py_func_dict_cache[key]
        except KeyError:
            # We don't have a cached version, so we need to generate it.
            py_func_dict = {}
            self._py_func_dict_cache[key] = py_func_dict
            for func_name, body in self._dynamic_funcs_python.items():
                if body != None:
                    exec body in self.namespace, locals()
                    py_func_dict[func_name] = locals()[func_name]

        # Add all the functions to our Network.
        for func_name, func in py_func_dict.items():
            setattr(self, func_name, func)
            self.namespace[func_name] = func

        if disable_c:
            return

        if curr_c_code is None:
            curr_c_code = self.get_c_code()
        # Search the cache.
        try:
            c_module = self._c_module_cache[curr_c_code]
        except KeyError:
            # Regenerate if needed.
            # Write C to file.
            module_name = self.output_c(curr_c_code)
            try:
                # Run f2py on the C. This may raise an exception if the command
                # fails.
                self.run_f2py(module_name, hide_f2py_output=True)
                c_module = __import__(module_name)
                if del_c_files:
                    os.unlink('%s.pyf' % module_name)
                    os.unlink('%s.c' % module_name)
                    try:
                        os.unlink('%s.so' % module_name)
                    except OSError:
                        pass
                    try:
                        os.unlink('%s.pyd' % module_name)
                    except OSError:
                        pass
            except ImportError, X:
                # Compiling C failed.
                logger.warn('Failed to import dynamically compiled C module %s '
                            'for network %s.' % (module_name, self.get_id()))
                logger.warn(X)
                # We stored None for the c_module, so we don't repeatedly
                # try compiling the same bad C code.
                c_module = None
            self._c_module_cache[curr_c_code] = c_module

        # Now we add all the appropriate functions to our Network.
        if c_module is not None:
            self.import_c_funcs_from_module(c_module)

    def get_c_code(self):
        # Combine all our C functions into one block of code
        c_code = []
        c_code.append('#include <math.h>')
        c_code.append('#include <stdio.h>')
        c_code.append('#include <float.h>')
        c_code.append('#include "mtrand.h"')
        c_code.append('#define exponentiale M_E')
        c_code.append('#define pi M_PI')
        c_code.append('double max(double a, double b){')
        c_code.append('return a > b ? a : b;}')
        c_code.append('double min(double a, double b){')
        c_code.append('return a < b ? a : b;}')
        
        # Function prototypes
        for func_name, proto in self._prototypes_c.items():
            c_code.append(proto)
            c_code.append('')
        # Functions necessary for SBML math support
        for function, body in self._common_func_strs_c.items():
            c_code.append(body)
            c_code.append('')
        # Function definitions
        for func_name, body in self._func_defs_c.items():
            c_code.append(body)
            c_code.append('')
        # The dynamic functions for this network
        for func_name, body in self._dynamic_funcs_c.items():
            c_code.append(body)
            c_code.append('')

        c_code = '\n'.join(code for code in c_code if code != None)
        return c_code

    def output_c(self, c_code, mod_name=None):
        logger.debug('Outputting C for network %s.' % self.get_id())
        if mod_name is None:
            semi_unique = str(time.time())
            mod_name = '%s_%i_%s_%s_%s' % (self.get_id(), SloppyCell.my_rank,
                                           SloppyCell.my_host,semi_unique[::-1],
                                           os.getpid())
            mod_name = mod_name.replace('-', '_')
            # f2py fails if base filename has a '.' in it.
            mod_name = mod_name.replace('.', '_')


        # Write the C code to a file.
        c_fd = open('%s.c' % mod_name, 'w')
        c_fd.write(c_code)
        c_fd.close()

        # Read in the signatures that we'll fill in
        # use different sig file if derivs are disabled
        if self.deriv_funcs_enabled != True:
            pyf_base_filename = os.path.join(SloppyCell.__path__[0], 
                                             'ReactionNetworks',
                                             'f2py_signatures_no_derivs.pyf')
        else:
            pyf_base_filename = os.path.join(SloppyCell.__path__[0], 
                                             'ReactionNetworks',
                                             'f2py_signatures.pyf')

        pyf_base_fd = file(pyf_base_filename, 'r')
        pyf_base = pyf_base_fd.read()
        
        pyf_base_fd.close()

        N_dyn = len(self.dynamicVars)
        N_const = len(self.constantVars)
        N_alg = len(self.algebraicVars)
        
        # Fill in the signatures.
        # We begin by generating the signatures for the dres_dparam code.
        # Since the number of optimizable parameters varies in each network,
        # we cannot include these functions in the f2p_signatres template.
        dres_dparams_code = ''
        if self.deriv_funcs_enabled == True:
            for wrt_ii, wrt in enumerate(self.optimizableVars.keys()):
                func_name = 'dres_d' + wrt
                dres_dparams_code += '    subroutine ' + func_name + '(time, dynamicVars, yprime, constants, pd)\n'
                dres_dparams_code += '        double precision intent(in) :: time\n'
                dres_dparams_code += '        double precision intent(in), dimension(' + str(N_dyn) + ') :: dynamicVars\n'
                dres_dparams_code += '        double precision intent(in), dimension(' + str(N_dyn) + ') :: yprime\n'
                dres_dparams_code += '        double precision intent(in), dimension(' + str(N_const) + ') :: constants\n'
                dres_dparams_code += '        double precision intent(out), dimension(' + str(N_dyn) + ') :: pd\n'
                dres_dparams_code += '    end subroutine ' + func_name + '\n'

        # Now update the template with the code we just generated and the other
        # names and dimensions
        pyf_code = pyf_base % {'dres_dparams': dres_dparams_code,
                               'mod_name': mod_name,
                               'N_dyn': N_dyn,
                               'N_const': N_const,
                               'N_alg': N_alg,
                               'N_rt': self.len_root_func}
        
        # Write out the signature file.
        pyf_fd = open('%s.pyf' % mod_name, 'w')
        pyf_fd.write(pyf_code)
        pyf_fd.close()

        return mod_name

    def run_f2py(self, mod_name, hide_f2py_output=True):
        logger.debug('Running f2py for network %s.' % self.get_id())
        from numpy.distutils.exec_command import exec_command
        # Run f2py. The 'try... finally...' structure ensures that we stop
        #  redirecting output if there's an exection in running f2py
        # These options assume we're working with mingw.
        win_options = ''
        if sys.platform == 'win32':
            win_options = '--compiler=mingw32 --fcompiler=gnu'
        try:
            if hide_f2py_output:
                redir = Utility.Redirector_mod.hideStdout()
                redir.start()
            sc_path = os.path.join(SloppyCell.__path__[0], 'ReactionNetworks')
            command = '-c %(win_options)s %(mod_name)s.pyf '\
                    '%(mod_name)s.c %(sc_path)s/mtrand.c -I%(sc_path)s'\
                    % {'win_options': win_options, 'mod_name': mod_name, 
                       'sc_path': sc_path}

            prev_cflags = os.getenv('CFLAGS', '')
            os.putenv('CFLAGS', prev_cflags + ' -Wno-unused-variable')
            # f2py wants an extra argument at the front here. It's not actually
            # used though...
            oldargv = sys.argv
            sys.argv =  ['f2py'] + command.split()
            import numpy.f2py.f2py2e
            output = numpy.f2py.f2py2e.run_compile()
            sys.argv = sys.argv
            os.putenv('CFLAGS', prev_cflags)
        except SystemExit, X:
            logger.warn('Call to f2py failed for network %s.' % self.get_id())
            logger.warn(X)
        finally:
            if hide_f2py_output:
                redir.stop()

    def import_c_funcs_from_module(self, module):
        for function in self._dynamic_funcs_c.keys():
            setattr(self, function, getattr(module, function))

    def takeDerivative(self, input, wrt, vars_used=None, simplify=True):
        """
        Take the derivative of a math expression wrt a given variable id.

        Does the chain rule through assigned variables.
        """
        output = ExprManip.diff_expr(input, wrt)

        if vars_used is None:
            vars_used = ExprManip.extract_vars(input)

        # What other assigned variables does input depend on?
        assigned_used = vars_used.difference(sets.Set([wrt]))
        assigned_used.intersection_update(sets.Set(self.assignedVars.keys()))
        # Do the chain rule for those variables
        for id in assigned_used:
            rule = self.assignmentRules.getByKey(id)
            d2 = self.takeDerivative(rule, wrt, simplify=False)
            if d2 != '0':
                d = ExprManip.diff_expr(input, id)
                output += ' + (%s) *(%s)' % (d, d2)

        # What other constant variables does input depend on?
        constant_used = vars_used.difference(sets.Set([wrt]))
        constant_used.intersection_update(sets.Set(self.constantVars.keys()))
        # Do the chain rule for those variables
        for id in constant_used:
            ic = self.get_var_ic(id)
            if isinstance(ic, str):
                d2 = self.takeDerivative(ic, wrt, simplify=False)
                if d2 != '0':
                    d = ExprManip.diff_expr(input, id)
                    output += ' + (%s) *(%s)' % (d, d2)

        if simplify and output != '0':
            return ExprManip.simplify_expr(output)
        else:
            return output

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

        # We can't copy the functions themselves. So we strip them out.
        for func in self._dynamic_funcs_python.keys():
            odict[func] = None
        for func in self._dynamic_funcs_c.keys():
            odict[func] = None
        odict['namespace'] = None
        # Let's not pickle these since they can be large and it would slow
        #  down parallel execution.
        odict['trajectory'] = None
        odict['ddv_dpTrajectory'] = None

        return odict

    def __setstate__(self, newdict):
        self.__dict__.update(newdict)
        self.namespace = copy.copy(self._common_namespace)

        # Recreate our namespace
        for func_id, func_str in self._func_strs.items():
            self.namespace[func_id] = eval(func_str, self.namespace, {})

        if not self._manualCrossReferences_flag:
            self._makeCrossReferences()
        if self.compiled:
            self.exec_dynamic_functions(self._last_disabled_c)

    def get_component_name(self, id, TeX_form=False):
        """
        Return a components's name if it exists, else just return its id.
        """
        # These are all the things that have names (except the network itself)
        complists = [self.variables, self.reactions, self.functionDefinitions,
                     self.events, self.constraints]
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
        # This was used to interface with PyDSTool.
        out = {}
        out['odes'] = dict(self.diff_eq_rhs.items())
        out['functions'] = {}
        for func_id, func_def in self.functionDefinitions.items():
            vars = ', '.join(func_def.variables)
            out['functions']['%s(%s)' % (func_id, vars)] = func_def.math
        out['parameters'] = dict([(id, var.value) for (id, var) 
                                  in self.constantVars.items()])
        out['assignments'] = dict(self.assignmentRules.items())
        out['events'] = dict([(event.trigger, 
                               dict(event.event_assignments.items()))
                              for event in self.events])
        out['constraints'] = dict([constraint.trigger for constraint in self.constraints])

        return out


    # Deprecated functions below.
    def addCompartment(self, id, size=1.0, name='', 
                       typicalValue=False,isConstant=True, isOptimizable=False):
        self.add_compartment(id = id, initial_size = size, name = name,
                             typical_value = typicalValue, 
                             is_constant = isConstant, 
                             is_optimizable = isOptimizable)

    addVariable = _add_variable
    
    def addSpecies(self, id, compartment, initialConcentration=None,
                   name='', typicalValue=None, is_boundary_condition=False,
                   isConstant=False, isOptimizable=False, uniprot_ids=None):
        self.add_species(id = id, compartment = compartment, 
                         initial_conc = initialConcentration, 
                         name = name,
                         typical_value = typicalValue, 
                         is_boundary_condition = is_boundary_condition, 
                         is_constant = isConstant,
                         is_optimizable = isOptimizable,
                         uniprot_ids=uniprot_ids)

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

    def addConstraint(self, id, trigger, message=None, name=''):
        self.add_constraint(id = id, trigger = trigger, message = message,
                            name = name)


    addFunctionDefinition = add_func_def
    addAssignmentRule = add_assignment_rule



def _exec_dynamic_func(obj, func, in_namespace={}, bind=True):
    """
    Create the executable function corresponding to func's functionBody.

    This only exists now for Trajectory_mod. It's not used by the Network
    object.
    """
    try:
        function_body = obj._dynamic_funcs_python.get(func)
    except (KeyError, AttributeError):
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
