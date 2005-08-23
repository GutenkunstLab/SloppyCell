import os
import copy
import types

import scipy

from SloppyCell import _TEMP_DIR, KeyedList
import Parsing

class Trajectory:
    def __init__(self, net, keyToColumn = None):
        if keyToColumn is not None:
            self.keyToColumn = keyToColumn
        else:
            keys = net.dynamicVars.keys + net.assignedVars.keys()
            self.keyToColumn = KeyedList(zip(keys, range(len(keys))))

        self.timepoints = scipy.zeros(0, scipy.Float)
	self.values = scipy.zeros((0, len(keyToColumn)), scipy.Float)

        self.assignedVarKeys = net.assignedVars.keys()
        self.dynamicVarKeys = net.dynamicVars.keys()
        self.optimizableVarKeys = net.optimizableVars.keys()

        self.constantVarValues = KeyedList([(id, var.value) for (id, var) in
                                            net.constantVars.items()])
        self.typicalVarValues = KeyedList([(id, var.typicalValue) for (id, var)
                                           in net.variables.items()])

        self.make_DoAssignmentInRow(net)
        if len(self.dynamicVarKeys + self.assignedVarKeys) < len(keyToColumn):
            self.make_DoSensAssignmentInRow(net)

    def append(self, other):
        if self.keyToColumn != other.keyToColumn:
            raise ValueError, 'Trajectories in append have different column keys!'
        if self.constantVarValues != other.constantVarValues:
            print 'WARNING: Constant variable values differ between appended trajectories!'

        if self.timepoints[-1] > other.timepoints[0]:
            print 'WARNING: Appending trajectory with earlier timepoints!'

        self.timepoints = scipy.concatenate((self.timepoints, other.timepoints))
        self.values = scipy.concatenate((self.values, other.values))

    def getVariableTrajectory(self, id):
        if not self.keyToColumn.has_key(id):
            raise ValueError, 'Variable %s not found in trajectory. Is it a constant variable?'
        return self.values[:, self.keyToColumn.getByKey(id)]

    def last_dynamic_var_values(self):
        """
        Return a list of the dynamic variable values at the last timepoint in
        the trajectory.
        """
        return [values[-1, self.keyToColumn[dv_id]] for dv_id in 
                self.dynamicVarKeys]

    def make_DoAssignmentInRow(self, net):
        functionBody = 'def DoAssignmentInRow(self, row, time):\n\t'

        if len(net.assignmentRules) > 0:
            for id, rule in net.assignmentRules.items():
                lhs = self.substituteVariableNames(id)
                rhs = net.substituteFunctionDefinitions(rule)
                rhs = self.substituteVariableNames(rhs)
                functionBody += lhs + ' = ' + rhs + '\n\t'

            functionBody = functionBody[:-2]
        else:
            functionBody += 'pass'

        f = file(os.path.join(_TEMP_DIR, 'DoAssignmentInRow.py'), 'w')
        print >> f, functionBody
        f.close()
        self.DoAssignmentInRow_functionBody = functionBody
        self.DoAssignmentInRow = None

    def make_DoSensAssignmentInRow(self, net):
        functionBody = 'def DoSensAssignmentInRow(self, row, time):\n\t'

        if len(net.assignmentRules) > 0:
	    for id, rule in net.assignmentRules.items():
                rule = net.substituteFunctionDefinitions(rule)
                derivWRTdv = {}
                for wrtId in net.dynamicVars.keys():
                    deriv = net.takeDerivative(rule, wrtId)
                    if deriv != '0':
                        derivWRTdv[wrtId] = deriv

		for optId in net.optimizableVars.keys():
                    lhs = self.substituteVariableNames('%s__derivWRT__%s' %
                                                           (id,optId))
                    rhs = '0'
                    # get derivative of assigned variable w.r.t.
                    #  dynamic variables
                    for wrtId, deriv in derivWRTdv.items():
                        rhs += ' + (%s) * %s__derivWRT__%s' % \
                                (deriv, wrtId, optId)

		    # now partial derivative w.r.t. optId
		    derivWRTp = net.takeDerivative(rule, optId)
		    if derivWRTp != '0' :
			rhs += ' + %s' % derivWRTp

                    if rhs != '0':
                        rhs = self.substituteVariableNames(rhs)
                        functionBody += '%s = %s\n\t' % (lhs, rhs)

            functionBody = functionBody[:-2]
        else:
            functionBody += 'pass'

        f = file(os.path.join(_TEMP_DIR, 'DoSensAssignmentInRow.py'), 'w')
        print >> f, functionBody
        f.close()
        self.DoSensAssignmentInRow_functionBody = functionBody
        self.DoSensAssignmentInRow = None

    def appendFromODEINT(self, timepoints, odeint_array):
        if self.DoAssignmentInRow is None:
            exec self.DoAssignmentInRow_functionBody
            self.DoAssignmentInRow = types.MethodType(DoAssignmentInRow, 
                                                      self, Trajectory)

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.keyToColumn)), 
                                  scipy.Float)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

	nDv = len(self.dynamicVarKeys)

        for ii, id in enumerate(self.keyToColumn.keys()[:nDv]):
            self.values[-numAdded:, self.keyToColumn.getByKey(id)] =\
                    odeint_array[:, ii]

        for time, row in zip(self.timepoints, self.values)[-numAdded:]:
            self.DoAssignmentInRow(row, time)

    def appendSensFromODEINT(self, timepoints, odeint_array):
        if self.DoAssignmentInRow is None:
            exec self.DoAssignmentInRow_functionBody
            self.DoAssignmentInRow = types.MethodType(DoAssignmentInRow, 
                                                      self, Trajectory)
        if self.DoSensAssignmentInRow is None:
            exec self.DoSensAssignmentInRow_functionBody
            self.DoSensAssignmentInRow = types.MethodType(DoSensAssignmentInRow,
                                                          self, Trajectory)

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.keyToColumn)), 
                                  scipy.Float)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

	nDv = len(self.dynamicVarKeys)

	# fill in trajectory
	for ii, id in enumerate(self.keyToColumn.keys()[0:nDv]):
            self.values[-numAdded:, self.keyToColumn.getByKey(id)] = \
                    odeint_array[:, ii]
	# ... and sensitivities
        for ii, dvId in enumerate(self.dynamicVarKeys):
            for jj, ovId in enumerate(self.optimizableVarKeys):
                self.values[-numAdded:, 
                            self.keyToColumn.getByKey((dvId, ovId))]\
                        = odeint_array[:, ii + (jj+1)*nDv]

        for time, row in zip(self.timepoints, self.values)[-numAdded:]:
            self.DoAssignmentInRow(row, time)
            self.DoSensAssignmentInRow(row, time)

    def __getstate__(self):
        odict = copy.copy(self.__dict__)
        odict['DoAssignmentInRow'] = None
        odict['DoSensAssignmentInRow'] = None
        return odict

    def substituteVariableNames(self, input):
	for id in Parsing.extractVariablesFromString(input):
            # convert it back to something keyToColumn will recognize
	    # had to use a form  dynVarName__derivWRT__optParamName for the
	    # sensitivity variable because otherwise,
            # Parsing.extractVariablesFromString gets confused
            splitId = id.split('__derivWRT__')
            if len(splitId) == 1:
	    	idname = splitId[0]
	    elif len(splitId) == 2:
	    	idname = tuple(splitId)
            else:
                raise 'Problem with id %s in Trajectory.substituteVariableNames'

	    if idname in self.keyToColumn.keys() :
                mapping = 'row[%i]' % self.keyToColumn.getByKey(idname)
            elif idname in self.constantVarValues.keys():
                # It must be a constantVars variable => substitute the
                #  value alone
                mapping = str(self.constantVarValues.getByKey(id))
            elif idname == 'time':
                mapping = 'time'
            else:
                raise 'Problem with idname %s in Trajectory.substituteVariableNames'

            input = Parsing.substituteVariableNamesInString(input, id, mapping)

        return input
