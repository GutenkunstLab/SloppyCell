"""
Created 6/2/2017

Author @Keeyan

Runs appropriate functions to test input data and model (XML edition
"""

import Utility
import time
import os
import sys
from pylab import *
import scipy
from SloppyCell.ReactionNetworks import *
from xml.etree import ElementTree as ET
import logging

logger = logging.getLogger('TestConstruct_XML')
logging.basicConfig()
this_module = sys.modules[__name__]
# Globally defined variables.  Should eventually be done away with.
step_factor = 300

# Here we define all of the modification functions for the networks
# ------------------------------------------------------------------------


def set_constant(network, action_root):
    for child in action_root:
        attributes = child.attrib
        network.set_var_constant(**attributes)
    return network


def add_species(network, action_root):
    for child in action_root:
        attributes = child.attrib
        network.add_species(**attributes)
    return network


def add_assignment(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        var_id = attributes.pop('id')
        network.add_assignment_rule(var_id, **attributes)
    return network


def set_optimizable(network, action_root):
    for child in action_root:
        attributes = child.attrib
        network.set_var_optimizable(**attributes)
    return network


def set_initial(network, action_root):
    for child in action_root:
        attributes = child.attrib
        network.set_initial_var_value(**attributes)
    return network


def start_fixed(network, action_root):
    attributes = action_root.attrib
    try:
        min = int(attributes['min'])
    except KeyError:
        min = 0
    try:
        max = int(attributes['max'])
    except KeyError:
        max = 1000
    Dynamics.integrate(network, [min, max])
    fp = Dynamics.dyn_var_fixed_point(network)
    network.set_dyn_var_ics(fp)
    return network


def add_event(network, action_root):
    # TODO: This could be easier to use
    for child in action_root:
        attributes = child.attrib
        assignment_dict = {}
        for assignment in child:
            attribs = assignment.attrib
            assignment_dict[attribs['id']] = attribs['func']
        network.add_event(event_assignments=assignment_dict, **attributes)
    return network


def add_parameter(network, action_root):
    for child in action_root:
        attributes = child.attrib
        network.add_parameter(**attributes)
    return network


def add_rate_rule(network, action_root):
    for child in action_root:
        attributes = child.attrib
        var_id = attributes.pop('id')
        network.add_rate_rule(var_id, **attributes)
    return network


def add_reaction(network, action_root):
    for child in action_root:
        pass_list = []
        attributes = child.attrib
        attributes = attributes.copy()
        var_id = attributes.pop('id')
        if 'stoich_var' in attributes and 'stoich_val' in attributes:
            stoich = {}
            try:
                value_s = int(attributes['stoich_val'])
            except ValueError:
                value_s = attributes['stoich_val']
            stoich[attributes['stoich_var']] = value_s
            pass_list.apppend(stoich)
        try:
            reaction_name = attributes.pop('reaction')
            try:
                kinetic_law = getattr(globals()['Reactions'],reaction_name)
                network.addReaction(kinetic_law, var_id, *pass_list, **attributes)
                return network
            except Exception as e:
                print e
                logger.warn("Please check that the reaction name is typed correctly")
        except KeyError:
            # No reaction specified
            kinetic_law = attributes.pop('kinetic_law')
            network.addReaction(var_id, kinetic_law, *pass_list, **attributes)
            return network
        for key in attributes:
            pass_list.append(attributes[key])
        network.addReaction(var_id, *pass_list)
    return network


def chop(word, character):
    if isinstance(character, str):
        if word.lower().find(character) > 0:
            return word.replace(character, '')
        else:
            return word
    else:
        if isinstance(character, list):
            for chara in character:
                if word.lower().find(chara) > 0:
                    return word.replace(chara, '')
                else:
                    pass
            return word


def routine_dict_drier(routine_dict):
    """
    Converts strings into the desired data type, allowing for easier use by other functions.
    
    :param routine_dict: A dictionary containing all the attributes from some XML element
    :return: Returns a copy of the original dictionary.  Directly modifying the original dictionary
             causes issues when writing back to file.
    """
    routine_dict_copy = routine_dict.copy()
    for key_r in routine_dict:
        parameter = routine_dict_copy[key_r]
        if isinstance(parameter, str):
            parameter_check = parameter.lower()
            # We check if the attribute is supposed to be a boolean once here and convert so we don't have to do it
            # in every subsequent function
            if parameter_check == "true":
                parameter = True
            elif parameter_check == "false":
                parameter = False
            else:
                try:
                    # Integer conversion will throw an error if there is a decimal
                    parameter = int(parameter)
                except ValueError:
                    # Check for possibility of a floating point number
                    try:
                        parameter = float(parameter)
                    except ValueError:
                        # The parameter is meant to be a string
                        # We allow certain attributes to be input as a comma-separated list.
                        if ',' in parameter:
                            placeholder = []
                            alpha = parameter.split(',')
                            for beta in alpha:
                                placeholder.append(beta.strip())
                            parameter = placeholder
        routine_dict_copy[key_r] = parameter
    return routine_dict_copy


def prior_drier(root):
    """
    Finds the appropriate values to put into a PriorInLog function to constrain the function.  
    Simple algebra tells us that 'val_p' is equal to the square root of the upper bound multiplied
    by the lower bound, and that 'x' is equal to the square root of the upper bound divided by 
    'val_p', or the square root of the 'val_p' divided by the lower bound.  
    The log of 'val_p' and 'x' become the parameters necessary for the PriorInLog function to bound
    the parameters correctly (the p value and sigma p value, respectively)

    :param root: Full dictionary output from ScellParser, contains the prior dictionary
    :return: A dictionary mapping each parameter to its respective (val_p,x) pair.  The (val_p,x) pair is 
      contained in a dictionary for convenience 
    """
    prior_dictionary = {}
    try:
        for child in root.iter('parameter'):
            for p in child.iter('prior'):
                bounds = p.attrib
                parameter = child.attrib.get('id')
                lower = float(bounds['lower'])
                upper = float(bounds['upper'])
                val_p = sqrt(lower * upper)
                x = sqrt(upper / val_p)
                # Just in case something goes wrong but doesn't cause any blatant errors
                assert (x == sqrt(val_p / lower))
                prior_dictionary[parameter] = {'val': val_p, 'x': x}
    except Exception as a:
        logger.warn('Prior Calculation malfunction: ensure prior bounds are numbers, and do not include things like'
                    '"log(x)", "sin(x)", "sqrt()"')
        logger.warn(a)
    return prior_dictionary


def add_residuals(root, model):
    prior_dictionary = prior_drier(root)
    for param in prior_dictionary.keys():
        values = prior_dictionary[param]
        res = Residuals.PriorInLog(param + '_prior', param, log(values['val']), log(values['x']))
        model.AddResidual(res)


def set_ic(fit_root, net):
    for child in fit_root.iter("parameter"):
        for ic in child.iter("ic"):
            original = child.attrib["id"]
            not_original = ic.attrib["id"]
            net.set_var_ic(not_original, original)


def find_vars(root, tag, key, value):
    tag_dict = dict()
    for child in root.iter(tag):
        tag_dict[child.attrib[key]] = child.attrib[value]
    return tag_dict


def construct_ensemble(params, m, steps=750, only_pruned=True, prune=0):
    j = m.jacobian_log_params_sens(log(params))  # use JtJ approximation to Hessian
    jtj = dot(transpose(j), j)
    print 'Beginning ensemble calculation.'
    ens, gs, r = Ensembles.ensemble_log_params(m, asarray(params), jtj, steps=steps)
    print 'Finished ensemble calculation.'
    if prune == 0:
        prune = int(math.ceil(steps/float(250)))
    pruned_ens = asarray(ens[::int(prune)])
    if only_pruned:
        return pruned_ens
    else:
        return ens, pruned_ens

        # print pruned_ens[:, 1]  # finds all values related to r3, goes by index


def plot_histograms(pruned_ens, keyed_list, ids='all', bins=10, log=True):
    if ids == 'all':
        ids = keyed_list.keys()
    if isinstance(ids, list):
        for param in ids:
            if param in keyed_list.keys():
                print "Plotting histogram for %s" % param
                param_index = keyed_list.index_by_key(param)
                fig = figure()
                fig.suptitle(param, fontweight='bold')
                param_array = pruned_ens[:, param_index]
                if log:
                    hist(scipy.log(param_array), normed=True, bins=bins)
                else:
                    hist(param_array, normed=True, bins=bins)
    else:
        param_index = keyed_list.index_by_key(ids)
        fig = figure()
        fig.suptitle(ids, fontweight='bold')
        param_array = pruned_ens[:, param_index]
        if log:
            hist(scipy.log(param_array), normed=True, bins=bins)
        else:
            hist(param_array, normed=True, bins=bins)


def plot_variables(pruned_ens = None, net = None, ids=None, time_r=65, start=0, points=100,
                   make_figure=True, color='g', bounds=(.025, .975), **kwargs):
    # TODO: make plots show up based on file input.
    times = linspace(start, time_r, points)
    if isinstance(bounds, list) or isinstance(bounds, str):
        bound_list = []
        if 'lower' in bounds:
            bound_list.append(.025)
        if 'upper' in bounds:
            bound_list.append(.975)
        if 'middle' in bounds:
            bound_list.append(.5)
        bounds = tuple(bound_list)
        if len(bounds) == 0:
            bounds = (.025, .975)
    if pruned_ens is not None:
        traj_set = Ensembles.ensemble_trajs(net, times, pruned_ens)
        quantiles = Ensembles.traj_ensemble_quantiles(traj_set, bounds)
    else:
        traj_set = kwargs['traj_set']
        quantiles = Ensembles.traj_ensemble_quantiles(traj_set, bounds)
        #quantiles = kwargs['quantiles']

    if make_figure:
        figure()

    print "Plotting trajectory set for %s" % ids
    for quantile in quantiles:
        plot(times, quantile.get_var_traj(ids), color)

    return traj_set


def cost_lm(params, m, optimize=True, plot=False, iterations=20):
    """
    This is the optimization routine.  It runs a cost analysis on the parameters provided and outputs a graph
    if asked.  Will do nothing if optimize and plot are set to false.  This can only be the case if the optimized 
    parameters are loaded from file.
    
    :param params: The parameters to find the cost for
    :param m: The model the contains the parameters
    :param optimize: A flag that specifies whether the optimized parameters have been already loaded or not
    :param plot: Whether to plot the model results or not, if set to 'True' the cost MUST be calculated
    :param iterations: The max iterations of the optimization function
    :return: returns the optimized parameters
    """

    if optimize:
        initial_cost = m.cost(params)
        print 'Initial Cost:', initial_cost
        params = Optimization.fmin_lm_log_params(m, params, maxiter=int(iterations), disp=False)
        optimized_cost = m.cost(params)
        print 'Optimized cost:', optimized_cost
        print 'Optimized parameters:', params
    if plot:
        if not optimize:
            optimized_cost = m.cost(params)
            print 'Optimized cost:', optimized_cost
            print 'Optimized parameters:', params
        figure()
        Plotting.plot_model_results(m)
    return params


def time_extract(experiment):
    """
    Goes through the experiment provided and finds the times
    """
    time_array = []
    for key_t in experiment.data.keys():
        for species_key in experiment.data[key_t].keys():
            times = experiment.data[key_t][species_key].keys()
            # Experiments hold the time data as keys in a dictionary
            time_array.append(max(times))
    return int(max(time_array)) + 5  # Arbitrary number added, still unsure how x-axis is calculated for plots


def key_parameters(fit_root, tag="parameter"):
    """
    Makes a keyed_list of the parameters to be fit.
    
    :param fit_root: XML root of parameters that need to be fit
    :param tag: The tag of what's being searched for, in this case the parameter tag almost always
    :return: Returns a keyed_list to be used for ensemble construction
    """
    param_list = []
    for child in fit_root.iter(tag):
        param_list.append((child.attrib['id'], float(child.attrib['value'])))
    params = KeyedList(param_list)
    return params


def find_by_tag(root, tag, attr=True, return_node=False):
    # Mostly obsolete, but could be useful for streamlining xml file searches
    found = False
    for child in root:
        if child.tag == tag:
            found = True
            if attr:
                if return_node:
                    return found, child, child.attrib
                else:
                    return found, child.attrib
    return found, None, None


def save_to_temp(obj, file_name, xml_file, node, routine="temp_file"):
    """
    Saves a single object to a file and stores it in the temp folder.  
    
    :param obj: The object to be saved
    :param file_name: The name of the actual xml file
    :param xml_file:  The xml tree object, used to write changes back to the file
    :param node: The node which should be updated
    :param routine: The name of the action routine, 
    :return: Doesn't return anything.
    """
    print "Saving " + routine + " to file"
    filename = '..\\temp\\' + routine + "_" + str(int(math.ceil(time.time()))) + ".bp"
    Utility.save(obj, filename)
    node.set('path', filename)
    xml_file.write(file_name)


def check_to_save(routine_function):
    """
    Decorator for save/load functionality of action routines
    """

    # TODO: Identify changes to variables and automatically neglect to load from file
    # TODO: Allow multiple files to be chosen from, maybe given names and dynamically accessed from command line?

    def wrapper(root=None, routine=None, **kwargs):
        """
        Handles the common "path", "save", and "use_path" attributes of most action routines
        
        :param root: (XMLTree) The root of the XML tree
        :param routine: (String) Name of the action routine to be done
        :param kwargs: Extra arguments that will be passed through the decorator function mostly unused
        :return: the result of the associated routine function is returned.  If an object is loaded it is simply
                 passed down to the associated routine function.
        """
        try:
            routine_dict = root.find(routine).attrib
        except AttributeError:
            try:
                current_root = kwargs['current_root']
                routine_dict = current_root.attrib
            except AttributeError:
                routine_dict = root.attrib
        routine_dict = routine_dict_drier(routine_dict)
        use_file = False
        save_file = False

        try:
            # We pop to get rid of the attributes related to saving, but we could just as easily
            # add **kwargs to all the routine actions to deal with the extra arguments
            use_file = routine_dict.pop('use_path')
        except KeyError as X:
            logger.debug("Path use not specified")
            logger.debug(X)
        try:
            save_file = routine_dict.pop('save')
        except KeyError as X:
            logger.debug("Save attribute not specified, will only save if no path currently exists.")
            logger.debug(X)
        if 'path' in routine_dict:
            saved_file_path = routine_dict.pop('path')
            if os.path.isfile(saved_file_path):
                # A file path is specified
                if use_file:
                    # User wants to load the object
                    print 'Successfully loaded %s objects from file' % routine
                    loaded_object = Utility.load(saved_file_path)
                    # Now we call the routine function and pass in the loaded object
                    # We can just return it because we don't have to save anything
                    # and we don't need to do anything else in this decorator
                    return routine_function(loaded_object=loaded_object, routine_dict=routine_dict, **kwargs)
                else:
                    # the user specified a path, but did not want to use the file.
                    if save_file:
                        # The user wants to overwrite the current path, so we delete it
                        # Todo: Old files don't have to be deleted, and the option to do so can be included
                        try:
                            os.unlink(saved_file_path)
                            print "Replacing old %s file with new one" % routine
                        except OSError:
                            pass
            else:
                # No valid file path is specified, but the path attribute has been specified
                # In this case we save to file anyway, and overwrite the invalid file path
                # with a valid one.  No trouble because we don't delete and files in this case.
                save_file = True
                pass
        else:
            # Path attribute omitted
            # In this case we don't even try to load, or save.
            # If the user doesn't specify the path attribute, but specified to save, the path attribute will be
            # added anyway when the function attempts to save, so we don't
            # even have to flip any of the booleans.
            pass
        # We can now call the routine function, but we don't return it until we figure out if we want to save it.
        routine_object = routine_function(routine_dict=routine_dict, **kwargs)
        try:
            xml_file = kwargs['xml_file']
            file_name = kwargs['file_name']
        except KeyError:
            logger.warn("Cannot save object for %s" % routine)
            save_file = False
        if save_file:
            node = root.find(routine)
            if node is None:
                    node = current_root
            save_to_temp(routine_object, xml_file=xml_file, file_name=file_name, node=node, routine=routine)
        return routine_object

    return wrapper


@check_to_save
def optimization(routine_dict, model, params, **kwargs):
    """
    Action routine for optimization.  Checks if the decorator has loaded an object, and if not, runs the regular
    optimization procedure.
    :return: the optimized parameters are returned and passed back to the decorator
    """
    try:
        optimized_params = kwargs['loaded_object']
        return cost_lm(optimized_params, model, optimize=False, **routine_dict)
    except KeyError:
        return cost_lm(params, model, **routine_dict)


@check_to_save
def ensemble(routine_dict, result_dictionary, model, **kwargs):
    """
    Action routine for ensemble construction.  Checks if an ensemble has been loaded, and if not, checks if the 
    optimized parameters have been calculated.
    """
    try:
        pruned_ensemble = kwargs['loaded_object']
        return pruned_ensemble
    except KeyError:
        pass
    try:
        optimized_params = result_dictionary['optimization']
    except KeyError:
        # We can just pass the unoptimized parameters to the construction function in this case.
        # Eventually should just call its dependent.
        logger.warn("Optimized parameters not found for ensemble construction, using unoptimized parameters instead")
        optimized_params = kwargs['params']
    return construct_ensemble(optimized_params, model, **routine_dict)


def histogram_r(current_root, result_dictionary, params, **kwargs):
    """
    Handles the histogram graphing routine.  Multiple variables being graphed with the same attributes can be 
    placed in a comma-separated list, as opposed to making an individual tag for each.  
    
    :param current_root: Histograms variables are not contained in the 'histogram' tag, so we need the xml tree root
    :param result_dictionary: dictionary of returned values from other action routines
    :param params: a keyed_list of the parameters to be fit
    :param kwargs: any extra arguments that must be passed through are handled by this
    :return: Does not return anything.  Will be put into the return dictionary as 'None'  
    """
    # Shadows numpy 'histogram()' function without underscore
    # Histograms cannot be loaded from file, they're based entirely on the ensemble.
    # If you have the ensemble you can make the histogram pretty quickly
    # Graph objects could potentially be saved, but it's probably not worth it.
    try:
        pruned_ens = result_dictionary['ensemble']
    except KeyError:
        # We might just not have done ensemble construction yet.
        # TODO: Either organize actions before calling, or have actions call their dependents
        pass
    for child in current_root:
        for variable in child:
            # We don't run the save function so the attributes don't get dried in the decorator
            attributes = variable.attrib
            attributes = routine_dict_drier(attributes)

            plot_histograms(pruned_ens, params, **attributes)


@check_to_save
def ensemble_traj(current_root, result_dictionary, time_r, net, **kwargs):
    print "Graphing trajectories"
    try:
        traj_dict = kwargs['loaded_object']
        object_loaded = True
    except KeyError:
        object_loaded = False
        traj_dict = {}
    color_array = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    try:
        pruned_ens = result_dictionary['ensemble']
    except KeyError:
        # We haven't constructed an ensemble yet.
        # TODO: Same as histogram, need ensemble to be called before this in some way
        pass
    try:
        time_r = current_root.attrib['time']
    except KeyError:
        logger.debug("No time specified, using default")
        # User did not specify a time, so we use the estimated one from the experiment
    # Need to dig into the root to get the variable information


    for child_r in current_root:
        for child in child_r:
            graphing_set = []
            index = 0
            for variable in child:
                # We don't run the save function so the attributes don't get dried in the decorator
                # We do it here instead
                attributes = variable.attrib
                attributes = routine_dict_drier(attributes)
                if 'color' not in attributes.keys():
                    attributes['color'] = color_array[index]
                    index += 1
                # We have 8 colors to choose from
                # TODO: Different graphs should not continue through the list and instead reset back to the beginning
                if index > 7:
                    index = 0
                graphing_set.append(attributes)
            figure_bool = True
            for species_attributes in graphing_set:
                species_id = species_attributes['ids']
                if not object_loaded:
                    traj_set = plot_variables(pruned_ens, net, time_r=time_r, make_figure=figure_bool, **species_attributes)
                    # We can use this dictionary to recreate the graphs, and also allow future functions to easily
                    # access the trajectory sets of different species
                    traj_dict[species_id] = (traj_set, figure_bool)
                else:
                    traj_set = traj_dict[species_id][0]
                    figure_bool = traj_dict[species_id][1]
                    #quantiles = traj_dict[species_id][2]
                    plot_variables(traj_set=traj_set, make_figure=figure_bool, time_r=time_r,
                                   **species_attributes)

                if figure_bool:
                    figure_bool = False
    return traj_dict


def trajectory_integration(result_dictionary, **kwargs):
    pass


@check_to_save
def create_Network(current_root, routine_dict, sbml_reference, network_dictionary, network_func_dictionary, **kwargs):
    try:
        network = kwargs['loaded_object']
        print network.id
        return network
    except KeyError:
        pass
    network_name = routine_dict['id']
    if 'from_file' in routine_dict:
        if routine_dict['from_file']:
            network = IO.from_SBML_file(sbml_reference)
    elif 'copy' in routine_dict:
        try:
            network_to_copy = network_dictionary[routine_dict['copy']]
            network = network_to_copy.copy(network_name)
        except KeyError:
            # Should call recursively to see if we just haven't made the network meant to be copied yet.
            # Can do this by ensuring this function is not called on networks that already exist in the dictionary

            # This is the case that could cause network to be undefined
            logger.warn("No network named %s to copy" % routine_dict['copy'])
    else:
        # We want to create a completely new network.
        network = Network(network_name)
    # Now that we have a network, we should modify it according to the xml file
    for child in current_root:
        r_name = child.tag.lower()
        network = network_func_dictionary[r_name](network, child)
    return network


def make_happen(root, experiment, xml_file=None, file_name=None, sbml_reference = None):
    result_dictionary = dict()
    network_dictionary = dict()

    action_function_dictionary = {'optimization': optimization, 'ensemble': ensemble, 'histogram': histogram_r,
                           'ensembletrajectories': ensemble_traj}
    network_func_dictionary = {'set_constant': set_constant, 'add_species': add_species,
                               'add_assignment': add_assignment, 'set_optimizable': set_optimizable,
                               'set_initial': set_initial, 'start_fixed': start_fixed,
                               'add_event': add_event, 'add_parameter': add_parameter,
                               'add_rate_rule': add_rate_rule, 'add_reaction': add_reaction}
    # TODO: Clean this mess up
    # These actions should generally happen for any input #
    # They have few enough options that there is little customization to be done and they are essential to setting
    # up the model.
    net = IO.from_SBML_file(sbml_reference)  # Set model from SBML reference
    time_r = time_extract(experiment)
    model = Model([experiment], [net])  # Use experiment to create model
    try:
        fit_root = root.find("Parameters").find("Fit")
        set_ic(fit_root, net)
        params = key_parameters(fit_root)
    except AttributeError:
        pass

    add_residuals(root, model)

    ###
    action_root = root.find("Actions")
    network_root = root.find("Networks")
    if network_root is not None:
        for network in network_root:
            net_name = network.attrib['id']
            print net_name
            network_dictionary[net_name] = create_Network(root=network_root, routine=net_name, xml_file=xml_file,
                                                          file_name=file_name, current_root=network,
                                                          sbml_reference=sbml_reference,
                                                          network_dictionary=network_dictionary,
                                                          network_func_dictionary=network_func_dictionary)
    else:
        # No network element is established, so we should just use the default network
        pass
    print network_dictionary
    for child in action_root:
        try:
            r_name = child.tag.lower()
            result_dictionary[r_name] = action_function_dictionary[r_name](root=action_root, routine=child.tag,
                                                                    model=model, params=params, current_root=child,
                                                                    xml_file=xml_file, file_name=file_name, time_r=time_r,
                                                                    net=net, result_dictionary=result_dictionary)
        except KeyError:
            logger.warn("No function associated with action tag %s" % child.tag)
    show()

        #
        # plot_histograms(pruned, ['r3', 'tao'], params)
        # show()
