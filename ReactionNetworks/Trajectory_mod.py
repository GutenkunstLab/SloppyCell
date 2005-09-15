import os
import copy
import types

import scipy

# Expose function definitions that SBML wants to see.
log, log10 = scipy.log, scipy.log10
exp = scipy.exp
cos, sin, tan = scipy.cos, scipy.sin, scipy.tan
acos, asin, atan = scipy.arccos, scipy.arcsin, scipy.arctan
cosh, sinh, tanh = scipy.cosh, scipy.sinh, scipy.tanh
arccosh, arcsinh, arctanh = scipy.arccosh, scipy.arcsinh, scipy.arctanh
exponentiale, pi = scipy.e, scipy.pi

import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList
import Network_mod
import Parsing

class Trajectory:
    _known_structures = []
    _known_function_bodies = []
    _dynamic_funcs = ['_assignment', '_sens_assignment']

    def __init__(self, net, key_column=None, is_sens=False):
        if key_column is not None:
            self.key_column = key_column
        else:
            keys = net.dynamicVars.keys() + net.assignedVars.keys()
            if is_sens:
                keys += [(cname, pname) for cname in keys
                         for pname in net.optimizableVars.keys()]
            self.key_column = KeyedList(zip(keys, range(len(keys))))

        # These are the main storage
        self.timepoints = scipy.zeros(0, scipy.Float)
	self.values = scipy.zeros((0, len(self.key_column)), scipy.Float)

        self.var_keys = net.variables.keys()
        self.dynamicVarKeys = net.dynamicVars.keys()
        self.assignedVarKeys = net.assignedVars.keys()
        self.optimizableVarKeys = net.optimizableVars.keys()

        self.const_var_values = KeyedList([(id, var.value) for (id, var) in
                                            net.constantVars.items()])
        self.typical_var_values = KeyedList([(id, var.typicalValue)
                                             for (id, var)
                                             in net.variables.items()])

        # To avoid generating our function bodies every Trajectory creation, we
        #  keep a list of known structures in the class itself.
        curr_structure = (net._get_structure(), self.key_column)
        for ii, struct in enumerate(self._known_structures):
            if curr_structure == struct:
                (self._assignment_functionBody,
                 self._sens_assignment_functionBody) = \
                        self._known_function_bodies[ii]
                break
        else:
            # We didn't find our our structure, so we have to make the function
            #  bodies. We don't want to make the sensitivity body unless we have
            #  to.
            self._assignment_functionBody = self.make__assignment(net)
            if is_sens:
                self._sens_assignment_functionBody =\
                        self.make__sens_assignment(net)

            self._known_structures.append(copy.deepcopy(curr_structure))
            bodies = (self._assignment_functionBody, 
                      getattr(self, '_sens_assignment_functionBody', None))
            self._known_function_bodies.append(bodies)

    def keys(self):
        return self.var_keys

    def get_times(self):
        return self.timepoints

    def append(self, other):
        if self.key_column != other.key_column:
            raise ValueError, 'Trajectories in append have different column keys!'
        if self.const_var_values != other.const_var_values:
            print 'WARNING: Constant variable values differ between appended trajectories!'

        if self.timepoints[-1] > other.timepoints[0]:
            print 'WARNING: Appending trajectory with earlier timepoints!'

        self.timepoints = scipy.concatenate((self.timepoints, other.timepoints))
        self.values = scipy.concatenate((self.values, other.values))

    def get_var_traj(self, id):
        if self.key_column.has_key(id):
            return self.values[:, self.key_column.get(id)]
        elif self.const_var_values.has_key(id):
            return scipy.ones(len(self.timepoints), scipy.Float) *\
                    self.const_var_values.get(id)
        else:
            raise ValueError, 'Variable %s not found in trajectory.' % id

    def get_var_vals_index(self, index):
        out = KeyedList(zip(self.keys(), [None]*len(self.keys())))
        for key, col in self.key_column.items():
            out.set(key, self.values[index, col])
        out.update(self.const_var_values)

        return out

    def get_var_vals(self, time, epsilon = 1e-6):
        index = scipy.argmin(abs(self.timepoints - time))
        time_range = self.timepoints[-1] - self.timepoints[0]
        if abs(self.timepoints[index] - time)/time_range > epsilon:
            print 'Time %f requested, closest time stored in trajectory is %f.'\
                    % (time, self.timepoints[index])

        return self.get_var_vals_index(index)

    def make__assignment(self, net):
        functionBody = ['def _assignment(self, values, times, start, end):']

        if len(net.assignmentRules) > 0:
            for id, rule in net.assignmentRules.items():
                lhs = self._sub_var_names(id)
                rhs = net.substituteFunctionDefinitions(rule)
                rhs = self._sub_var_names(rhs)
                functionBody.append('%s = %s' % (lhs, rhs))
        else:
            functionBody.append('pass')

        return '\n\t'.join(functionBody) + '\n'

    def make__sens_assignment(self, net):
        functionBody = ['def _sens_assignment(self, values, times, start, end):'
                        ]

        if len(net.assignmentRules) > 0:
	    for id, rule in net.assignmentRules.items():
                rule = net.substituteFunctionDefinitions(rule)
                derivWRTdv = {}
                for wrtId in net.dynamicVars.keys():
                    deriv = net.takeDerivative(rule, wrtId)
                    if deriv != '0':
                        derivWRTdv[wrtId] = deriv

		for optId in net.optimizableVars.keys():
                    lhs = self._sub_var_names('%s__derivWRT__%s' % (id,optId))
                    rhs = []
                    # get derivative of assigned variable w.r.t.
                    #  dynamic variables
                    for wrtId, deriv in derivWRTdv.items():
                        rhs.append('(%s) * %s__derivWRT__%s' % 
                                   (deriv, wrtId, optId))

		    # now partial derivative w.r.t. optId
		    derivWRTp = net.takeDerivative(rule, optId)
		    if derivWRTp != '0':
			rhs.append(derivWRTp)

                    if rhs:
                        rhs = ' + '.join(rhs)
                        rhs = self._sub_var_names(rhs)
                        functionBody.append('%s = %s' % (lhs, rhs))
        else:
            functionBody.append('pass')

        return '\n\t'.join(functionBody) + '\n'

    def appendFromODEINT(self, timepoints, odeint_array):
        if getattr(self, '_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_assignment')

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.key_column)), 
                                  scipy.Float)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

        for ii, id in enumerate(self.dynamicVarKeys):
            self.values[-numAdded:, self.key_column.get(id)] =\
                    odeint_array[:, ii]

        self._assignment(self.values, self.timepoints, -numAdded, None)

    def appendSensFromODEINT(self, timepoints, odeint_array):
        if getattr(self, '_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_assignment')

        if getattr(self, '_sens_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_sens_assignment')

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.key_column)), 
                                  scipy.Float)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

	nDv = len(self.dynamicVarKeys)

	# fill in trajectory
	for ii, dvId in enumerate(self.dynamicVarKeys):
            self.values[-numAdded:, self.key_column.get(dvId)] = \
                    odeint_array[:, ii]
	# ... and sensitivities
        for ii, dvId in enumerate(self.dynamicVarKeys):
            for jj, ovId in enumerate(self.optimizableVarKeys):
                self.values[-numAdded:, 
                            self.key_column.get((dvId, ovId))]\
                        = odeint_array[:, ii + (jj+1)*nDv]

        self._assignment(self.values, self.timepoints, -numAdded, None)
        self._sens_assignment(self.values, self.timepoints, -numAdded, None)

    def __getstate__(self):
        odict = copy.copy(self.__dict__)
        odict['_assignment'] = None
        odict['_sens_assignment'] = None
        return odict

    def _sub_var_names(self, input):
	for id in Parsing.extractVariablesFromString(input):
            # convert it back to something key_column will recognize
	    # had to use a form  dynVarName__derivWRT__optParamName for the
	    # sensitivity variable because otherwise,
            # Parsing.extractVariablesFromString gets confused
            splitId = id.split('__derivWRT__')
            if len(splitId) == 1:
	    	idname = splitId[0]
	    elif len(splitId) == 2:
	    	idname = tuple(splitId)
            else:
                raise 'Problem with id %s in Trajectory._sub_var_names'

	    if idname in self.key_column.keys():
                mapping = 'values[start:end, %i]' % self.key_column.get(idname)
            elif idname in self.const_var_values.keys():
                # It must be a constantVars variable => substitute the
                #  value alone
                mapping = str(self.const_var_values.get(id))
            elif idname == 'time':
                mapping = 'times[start:end, %i]'
            else:
                raise 'Problem with idname %s in Trajectory._sub_var_names'

            input = Parsing.substituteVariableNamesInString(input, id, mapping)

        return input

    # Deprecated
    def last_dynamic_var_values(self):
        """
        Return a list of the dynamic variable values at the last timepoint in
        the trajectory.
        """
        return [self.values[-1, self.key_column.get(dv_id)] for dv_id in 
                self.dynamicVarKeys]

    getVariableTrajectory = get_var_traj
