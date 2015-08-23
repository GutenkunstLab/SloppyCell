from __future__ import division

import logging
logger = logging.getLogger('RxnNets.Trajectory_mod')

import os
import copy
import types

import scipy
import scipy.interpolate

import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList
import Network_mod
import SloppyCell.ExprManip as ExprManip

class Trajectory:
    _known_structures = []
    _known_function_bodies = []
    _dynamic_funcs = ['_assignment', '_sens_assignment']
    # These are-predefined functions we want in our working namespace
    _common_namespace = {'log': scipy.log,
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
                         'pow': scipy.power,
                         'sqrt': scipy.sqrt,
                         'exponentiale': scipy.e,
                         'pi': scipy.pi,
                         'min': scipy.minimum,
                         'max': scipy.maximum
                         }

    def __init__(self, net, key_column=None, is_sens=False, holds_dt=False,
                 empty=False, const_vals=None):
        if empty:
            return

        if key_column is not None:
            self.key_column = key_column
        else:
            keys = net.dynamicVars.keys() + net.assignedVars.keys()
            if is_sens:
                keys.extend([(cname, pname) for cname in keys
                             for pname in net.optimizableVars.keys()])
            if holds_dt:
                for keyname in copy.copy(keys):
                    if isinstance(keyname,str):
                        keys.append((keyname,'time'))
                    else: # key is a tuple
                        keys.append(keyname + ('time',))

            self.key_column = KeyedList(zip(keys, range(len(keys))))

        # These are the main storage
        self.timepoints = scipy.zeros(0, scipy.float_)
	self.values = scipy.zeros((0, len(self.key_column)), scipy.float_)

        self.var_keys = net.variables.keys()
        self.dynamicVarKeys = net.dynamicVars.keys()
        self.assignedVarKeys = net.assignedVars.keys()
        self.optimizableVarKeys = net.optimizableVars.keys()

        # We do an 'evaluate_expr' here to take care of constant variables that 
        #  are initialized by other variables
        if const_vals is None:
            self.const_var_values = KeyedList([(id, net.evaluate_expr(id)) for 
                                               id in net.constantVars.keys()])
        else:
            self.const_var_values = KeyedList(zip(net.constantVars.keys(), 
                                                  const_vals))

        self.typical_var_values = KeyedList([(id, var.typicalValue)
                                             for (id, var)
                                             in net.variables.items()])
        self.event_info = None
        self.tcks = {} # need this to store the interpolation information
        self.dytcks = {} # to store interpolated vector field info

        # We make a copy of the Network's namespace.
        self._func_strs = copy.copy(net._func_strs)
        self.namespace = copy.copy(self._common_namespace)
        for func_id, func_str in self._func_strs.items():
            self.namespace[func_id] = eval(func_str, self.namespace, {})

        # To avoid generating our function bodies every Trajectory creation, we
        #  keep a list of known structures in the class itself.
        curr_structure = (net._last_structure, self.key_column, is_sens)
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
            self._assignment_functionBody = self._make__assignment(net)
            if is_sens:
                self._sens_assignment_functionBody =\
                        self._make__sens_assignment(net)

            self._known_structures.append(copy.deepcopy(curr_structure))
            bodies = (self._assignment_functionBody, 
                      getattr(self, '_sens_assignment_functionBody', None))
            self._known_function_bodies.append(bodies)

    def __len__(self):
        return len(self.timepoints)

    def __getitem__(self, this_slice):
        # XXX: This is very memory ineffcient. A whole copy is made and then
        #      items discarded. Should really fix this.
        new_traj = copy.deepcopy(self)
        if not isinstance(this_slice, slice):
            # This is generally just the case where we have a single index.
            this_slice = slice(this_slice, this_slice+1, None)
        new_traj.timepoints = self.timepoints[this_slice]
        new_traj.values = self.values[this_slice]

        # XXX: For now just clear out the interoplating info. This could
        #      be handled cleanly in the future
        if self.tcks or new_traj.dytcks:
            logger.warn('Interpolating functions must be recreated after '
                        'slicing a trajectory. Could be fixed.')
            new_traj.tcks = {} 
            new_traj.dytcks = {} 

        return new_traj

    def time_slice(self, start, stop, eps=1e-6):
        """
        Return a new trajectory containing only times from start to stop.
        """
        start_index = self._get_time_index(start, eps)
        stop_index = self._get_time_index(stop, eps)
        return self[start_index:stop_index+1]

    def keys(self):
        return self.var_keys

    def get_times(self):
        return self.timepoints

    def add_event_info(self, net, eventinfo, time_traj_ended, include_extra_event_info = False):
        """
        Add information about the network state at each event execution
        """

        # if include_assigned_vals is True, then append a 4th item to
        # the eventinfo tuple that contains the values of the assigned
        # vars in the network, evaluated at each event state

        if include_extra_event_info == True:
            (te, ye_pre, ye_post, ie) = eventinfo
            assigned_states = []
            prev_vals = net.get_var_vals()
            for t, y, ii in zip(te, ye_pre, ie):
                net.updateVariablesFromDynamicVars(y, t)
                a_ids = [id for id in net.assignedVars.keys()]
                a_vals = [net.get_var_val(id) for id in net.assignedVars.keys()]
                # we store the assigned variables in a dictionary for easy retrieval later
                a_state = dict(zip(a_ids,a_vals))
                assigned_states.append(a_state)
            net.set_var_vals(prev_vals, time_traj_ended)
            eventinfo = (te,ye_pre,ye_post,ie,assigned_states)

        # the information is a triple (te,ye,ie), or a quintuple
        # if include_extra_event_info is True
        self.event_info = eventinfo  

        

    def append(self, other):
        if self.key_column != other.key_column:
            raise ValueError, 'Trajectories in append have different column keys!'
        if self.const_var_values != other.const_var_values:
            logger.warn('Constant variable values differ between appended trajectories!')

        if self.timepoints[-1] > other.timepoints[0]:
            logger.warn('Appending trajectory with earlier timepoints!')

        self.timepoints = scipy.concatenate((self.timepoints, other.timepoints))
        self.values = scipy.concatenate((self.values, other.values))

    def get_var_typical_val(self, id):
        return self.typical_var_values.get(id)

    def get_var_traj(self, id):
        if self.key_column.has_key(id):
            return self.values[:, self.key_column.get(id)]
        elif self.const_var_values.has_key(id):
            return scipy.ones(len(self.timepoints), scipy.float_) *\
                    self.const_var_values.get(id)
        elif id == 'time':
            return self.get_times()
        elif len(id) == 2 and id[0] in self.const_var_values.keys():
            # Requesting sensitivity of a constant variable.
            return 0*self.get_times()
        else:
            raise ValueError, 'Variable %s not found in trajectory.' % str(id)

    def _get_time_index(self, time, eps=1e-6):
        """
        Return the index of the stored value closest to the requested time.

        Prints a warning if the difference between the requested time and the
        stored time is greater than a fraction eps of the trajectory length.
        """
        index = scipy.argmin(abs(self.timepoints - time))
        time_range = self.timepoints[-1] - self.timepoints[0]
        if abs(self.timepoints[index] - time)/time_range > eps:
            logger.warn('Time %f requested, closest time stored in trajectory '
                        'is %f.' % (time, self.timepoints[index]))
        return index

    def get_var_vals(self, time, eps=1e-6):
        """
        Return a KeyedList of the values of the trajectory's variables at the
        given time.

        Prints a warning if the difference between the requested time and the
        stored time is greater than a fraction eps of the trajectory length.
        """
        index = self._get_time_index(time, eps)
        return self.get_var_vals_index(index)

    def get_var_val(self, var_id, time, eps=1e-6):
        """
        Return the value of the given variable at the given time.

        Prints a warning if the difference between the requested time and the
        stored time is greater than a fraction eps of the trajectory length.
        """
        index = self._get_time_index(time, eps)
        return self.get_var_val_index(var_id, index)

    def get_var_vals_index(self, index):
        """
        Return a KeyedList of the values of the trajectory's variables at the
        given index.
        """
        out = KeyedList([(key, self.get_var_val_index(key, index)) for
                          key in ['time'] + self.keys()])
        return out

    def get_var_val_index(self, var_id, index):
        """
        Return the value of the given variable at the given index.
        """
        if self.key_column.has_key(var_id):
            col = self.key_column.get(var_id)
            return self.values[index, col]
        elif self.const_var_values.has_key(var_id):
            return self.const_var_values.get(var_id)
        elif var_id == 'time':
            return self.timepoints[index]
        elif len(var_id) == 2 and var_id[1] in self.const_var_values.keys():
            # Requesting sensitivity of a constant variable
            return 0

    def get_dynvar_vals_index(self, index):
        """
        Return a KeyedList of the values of the trajectory's dynamic variables
        at the given index.
        """
        out = KeyedList([(key, self.get_var_val_index(key, index)) for
                          key in self.dynamicVarKeys])
        return out

    def get_dynvar_vals(self, time, eps=1e-6):
        """
        Return a KeyedList of the values of the trajectory's dynamic variables
        at the given time.

        Prints a warning if the difference between the requested time and the
        stored time is greater than a fraction eps of the trajectory length.
        """
        index = self._get_time_index(time, eps)
        return self.get_dynvar_vals_index(index)

    def _make__assignment(self, net):
        functionBody = ['def _assignment(self, values, times, start, end):']

        for id in self.const_var_values.keys():
            functionBody.append("%s = self.const_var_values.get('%s')" % 
                                (id, id))

        if len(net.assignmentRules) > 0:
            for id, rule in net.assignmentRules.items():
                lhs = self._sub_var_names(id)
                #rhs = net.substituteFunctionDefinitions(rule)
                rhs = self._sub_var_names(rule)
                functionBody.append('# Assignment rule %s = %s' % (id, rule))
                functionBody.append('%s = %s' % (lhs, rhs))
        else:
            functionBody.append('pass')

        return '\n\t'.join(functionBody) + '\n'

    def _make__sens_assignment(self, net):
        functionBody = ['def _sens_assignment(self, values, times, start, end):'
                        ]

        for id in self.const_var_values.keys():
            functionBody.append("%s = self.const_var_values.get('%s')" % 
                                (id, id))

        if len(net.assignmentRules) > 0:
	    for id, rule in net.assignmentRules.items():
                #rule = net.substituteFunctionDefinitions(rule)
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

    def appendFromODEINT(self, timepoints, odeint_array, holds_dt = False):
        if getattr(self, '_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_assignment',
                                           self.namespace, bind=False)

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.key_column)),
                                  scipy.float_)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

        for ii, id in enumerate(self.dynamicVarKeys):
            self.values[-numAdded:, self.key_column.get(id)] =\
                    odeint_array[:, ii]

        self._assignment(self.values, self.timepoints, -numAdded, None)
        if holds_dt :
            for ii, id in enumerate(self.dynamicVarKeys) :
                self.values[-numAdded:, self.key_column.get((id,'time'))] = \
                    odeint_array[:,ii+len(self.dynamicVarKeys)]

    def appendSensFromODEINT(self, timepoints, odeint_array, holds_dt = False):
        if getattr(self, '_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_assignment',
                                           self.namespace, bind=False)

        if getattr(self, '_sens_assignment', None) is None:
            Network_mod._exec_dynamic_func(self, '_sens_assignment',
                                           self.namespace, bind=False)

        numAdded = odeint_array.shape[0]
        addedValues = scipy.zeros((numAdded, len(self.key_column)),
                                  scipy.float_)

        self.values = scipy.concatenate((self.values, addedValues))
        self.timepoints = scipy.concatenate((self.timepoints, timepoints))

        nDv = len(self.dynamicVarKeys)
        nOv = len(self.optimizableVarKeys)

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

        if holds_dt :
        # fill in the time derivative of the trajectory
            for ii, dvId in enumerate(self.dynamicVarKeys):
                self.values[-numAdded:, self.key_column.get((dvId,'time'))] = \
                    odeint_array[:, ii + nDv*(nOv+1)]
        # ... and of the sensitivities
            for ii, dvId in enumerate(self.dynamicVarKeys):
                for jj, ovId in enumerate(self.optimizableVarKeys):
                    self.values[-numAdded:,
                            self.key_column.get((dvId, ovId,'time'))]\
                        = odeint_array[:, ii + (jj+1)*nDv + nDv*(nOv+1)]

    def copy_subset(self, keys=None):
        """
        Return a copy of this trajectory containing only the variables specified
        in keys.

        If keys is None, all variables are included.
        """
        if keys is None:
            keys = self.key_column.keys()

        state = self.__getstate__()

        # Only need those keys that are stored in the values array. The
        #  rest will be copied easily.
        keys = [key for key in keys if self.key_column.has_key(key)]
        new_key_column = KeyedList(zip(keys, range(len(keys))))
        state['key_column'] = new_key_column

        new_values = scipy.zeros((len(self.values), len(new_key_column)), 
                                 scipy.float_)
        for key, new_col in new_key_column.items():
            old_col = self.values[:, self.key_column.get(key)]
            new_values[:, new_col] = old_col.copy()
        state['values'] = new_values

        new_traj = Trajectory(None, empty=True)
        new_traj.__setstate__(state)
        return new_traj

    def __getstate__(self):
        odict = copy.copy(self.__dict__)
        odict['_assignment'] = None
        odict['_sens_assignment'] = None
        odict['namespace'] = None
        return odict

    def __setstate__(self, newdict):
        self.__dict__.update(newdict)
        # Remake our namespace
        self.namespace = copy.copy(self._common_namespace)
        for func_id, func_str in self._func_strs.items():
            self.namespace[func_id] = eval(func_str, self.namespace, {})

    def _sub_var_names(self, input):
        mapping_dict = {}
	for id in ExprManip.extract_vars(input):
            # convert it back to something key_column will recognize
	    # had to use a form  dynVarName__derivWRT__optParamName for the
	    # sensitivity variable because otherwise,
            # extract_vars gets confused
            splitId = id.split('__derivWRT__')
            if len(splitId) == 1:
	    	idname = splitId[0]
	    elif len(splitId) == 2:
	    	idname = tuple(splitId)
            else:
                raise 'Problem with id %s in Trajectory._sub_var_names' % id

	    if idname in self.key_column.keys():
                mapping = 'values[start:end, %i]' % self.key_column.get(idname)
            elif idname in self.const_var_values.keys():
                # Don't substitute for constant variable names. Those will
                #  be taken care of earlier in the method.
                continue
            elif idname == 'time':
                mapping = 'times[start:end]'
            else:
                raise 'Problem with idname %s in Trajectory._sub_var_names' % id
            mapping_dict[id] = mapping

        input = ExprManip.sub_for_vars(input, mapping_dict)

        return input

    def build_interpolated_traj(self) :
        """ Given that a trajectory exists, build_interpolated_traj will create 
        the coefficients for the spline interpolatation.
        The spline can then be evaluated using 
        Trajectory.evaluate_interpolated_traj or
        Trajectory.evaluate_interpolated_trajs """
        te,ye,ie = self.event_info[:3]
        teIndices = []

        if len(te) == 0 : # no events
            intervals = [(0,len(self.timepoints))]
        else :
            # At an event there are two time points in the trajectory that
            # are the same (=tevent) but we want the second one
            last_t = -1
            for tevent in te :
                if tevent != last_t:
                    teIndices.append(scipy.nonzero(self.timepoints==tevent)[0][1])
                    last_t = tevent                 

            # don't expect there to be an event at 0, if there is this will be
            # messed up
            teIndicesWith0 = list(teIndices[0:])
            teIndicesWith0.insert(0,0)
            # put in the last time point as well, again a problem if there's an
            # event at the last time
            teIndices.extend([len(self.timepoints)])
            intervals = zip(teIndicesWith0,teIndices)

        self.tcks = {}

        for (start_ind,end_ind) in intervals :
            start_time, end_time = self.timepoints[start_ind], \
                    self.timepoints[end_ind-1]
            curTimes = self.timepoints[start_ind:end_ind]
            k = min(5,end_ind-start_ind-1)
            ys = [self.get_var_traj(dv_id)[start_ind:end_ind]
                    for dv_id in self.key_column.keys()]

            self.tcks[(start_time,end_time)] = [scipy.interpolate.splrep(curTimes,scipy.asarray(y),k=k,s=0) for y in ys]

        #return self.tcks # do we want to return this?

    def evaluate_interpolated_trajs(self,time,subinterval,der=0) :
        """
        This is a version of evaluate_interpolated_traj that returns all the
        values of the dynamic variables and requires you to pass in the
        appropriate subinterval between events (that can be found in 
        Trajectory.tcks.keys() )
        Faster than calling evaluate_interpolated_traj repeatedly
        """
        local_tcks = self.tcks

        nDVs = len(self.dynamicVarKeys)
        dv_y = [scipy.interpolate.splev(time, local_tcks[subinterval][dv_ind],
                                        der=der) for dv_ind in range(0,nDVs)]

        return dv_y

    def evaluate_interpolated_traj(self,dv_id,time,subinterval=None,der=0) :
        """ Needs Trajectory.build_interpolated_traj() to be called first

        Arguments:
        dvid         the name of the component of the trajectory you wish to 
                     evaluate
        time         a vector of times or a scalar
        subinterval  an optional argument specifying the time interval 
                     between events that the time argument lies (but given a 
                     valid time, it will be found automatically)
        der          the derivative of the spline function you want, the order
                     of the derivative will be constrained by the order of the 
                     interpolated spline
        Outputs:
        A single scalar value (if time input is a scalar)
        or
        (returned_times, interpolated_trajectory at those times) if times is a
        vector

        Note: It is necessary to have a returned_times argument too, in case 
              the times passed in happens to have a timepoint that corresponds 
              to an event time, which often has two trajectory values associated
              with it.
        """
        if scipy.isscalar(time) :
            time = scipy.asarray([time]) # if a scalar was passed in, convert to an array
        else :
            time = scipy.asarray(time)
        local_tcks = self.tcks
        sorted_intervals = scipy.sort(local_tcks.keys(),axis=0)

        if subinterval is not None : # confine things to just one interval
            if subinterval not in local_tcks.keys() :
                raise "Not a valid subinterval (not in Trajectory.tcks.keys())"
            else :
                sorted_intervals = [[subinterval[0],subinterval[1]]]
                interval_start_ind = 0
                interval_end_ind = 0
        else :
            # sorted_intervals ends up being a list of lists, each length 2, not tuples anymore
            for interval_ind, interval in enumerate(sorted_intervals) :
                start_time, end_time = interval[0],interval[1]
                if (time[0] >= start_time) :
                    interval_start_ind = interval_ind
                if (time[-1] <= end_time) :
                    interval_end_ind = interval_ind
                    break

        dv_y = []
        returned_times = []
        dv_ind = self.key_column.keyToIndex[dv_id]
        for interval in sorted_intervals[interval_start_ind:(interval_end_ind+1)] :
            currTimes = scipy.compress( scipy.logical_and((time>=interval[0]),(time<=interval[1])) , time )
            startslice, endslice = 0, None
            if len(currTimes) > 1 :
                if (currTimes[0]==currTimes[1]) :
                # skip the first time point because it's repeated
                    startslice = 1
                if (currTimes[-1]==currTimes[-2]) :
                # skip the last time point because it's repeated
                    endslice = -1
                dv_y.extend( scipy.interpolate.splev(currTimes[startslice:endslice],
                                local_tcks[(interval[0],interval[1])][dv_ind],der=der) )
                returned_times.extend(currTimes[startslice:endslice])
            elif len(currTimes) == 1: # explicitly check, because len(currTimes) could = 0
                dv_y.extend( [ scipy.interpolate.splev(currTimes, local_tcks[(interval[0],interval[1])][dv_ind],der=der) ])
                returned_times.extend(currTimes[startslice:endslice])

        if len(returned_times) == 1 :
            return dv_y[0]
        else :
            return returned_times,dv_y

    def to_file(self, file_name, out_vars=None, separator=', '):
        """
        Output the given variables to a file.

        file_name   Name of the file to use for output (will be overwritten!)
        out_vars    List of variable ids ot output. If None, default is
                    'time', dynamic variables
        separator   The separator to use between columns.
        """
        if out_vars is None:
            out_vars = ['time'] + self.dynamicVarKeys
    
        first_line = separator.join(out_vars) + os.linesep
        f = file(file_name, 'w')
        f.write(first_line)
    
        out_array = []
        for var in out_vars:
            out_array.append(self.get_var_traj(var))
    
        out_array = scipy.transpose(out_array)
        scipy.savetxt(f, out_array, delimiter=separator)
        f.close()

    def merge(self, traj):
        """
        Merge one trajectory with another.
        """

        # Verify that both trajectories had the same keys in the same order
        for traj_key, traj_index in traj.key_column.items():
            if self.key_column.get(traj_key) != traj_index:
                raise ValueError('Trajectories are not mergeable')

        self.values = scipy.array(list(self.values)+list(traj.values))
        last_time = self.timepoints[-1]
        updated_times = traj.timepoints
        self.timepoints = scipy.array(list(self.timepoints)+list(updated_times))
        
        (te,ye,ie) = traj.event_info
        updated_event_times = te
        (self_te,self_ye,self_ie) = traj.event_info        

        self.event_info = (scipy.array(list(self_te)+list(updated_event_times)),
                           scipy.array(list(self_ye)+list(ye)),
                           scipy.array(list(self_ie)+list(ie)))

        self.const_var_values = traj.const_var_values



    # Deprecated
    def last_dynamic_var_values(self):
        """
        Return a list of the dynamic variable values at the last timepoint in
        the trajectory.
        """
        return [self.values[-1, self.key_column.get(dv_id)] for dv_id in
                self.dynamicVarKeys]

    getVariableTrajectory = get_var_traj

    try:
        import psyco
        psyco.bind(scipy.interpolate.splrep)
        psyco.bind(scipy.interpolate.splev)
        psyco.bind(evaluate_interpolated_trajs)
        psyco.bind(evaluate_interpolated_traj)
    except ImportError:
        pass
