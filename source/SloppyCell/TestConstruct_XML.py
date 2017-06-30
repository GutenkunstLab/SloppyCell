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


def merge_dicts(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def recursive_list_builder(node, order):
    if len(node) < 1:
        return order
    for child in node:
        order.append(child)
        return recursive_list_builder(child,order)

def hash_routine(root):
    root_iter = root.getiterator()
    hash_list = []
    for child in root_iter:
        hash_list.append((child.tag, (tuple(child.attrib.keys()), tuple(child.attrib.values()))))
        hash_list.sort()
    return hash(tuple(hash_list))


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

# Here we define all of the modification functions for the networks
# ------------------------------------------------------------------------


def set_constant(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.set_var_constant(**attributes)
    return network


def add_species(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.add_species(**attributes)
    return network


def add_assignment(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
        var_id = attributes.pop('id')
        network.add_assignment_rule(var_id, **attributes)
    return network


def set_optimizable(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.set_var_optimizable(**attributes)
    return network


def set_initial(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.set_initial_var_value(**attributes)
    return network


def start_fixed(network, action_root):
    attributes = action_root.attrib
    attributes = routine_dict_drier(attributes)
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
        attributes = routine_dict_drier(attributes)
        network.add_parameter(**attributes)
    return network


def add_rate_rule(network, action_root):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
        var_id = attributes.pop('id')
        network.add_rate_rule(var_id, **attributes)
    return network


def add_reaction(network, action_root):
    for child in action_root:
        pass_list = []
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
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


def prior_drier(root, model):
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
    all_params = model.get_params()
    all_params = all_params.items()

    try:
        for child in root.iter('parameter'):
            for p in child.iter('prior'):
                param_list = []
                bounds = p.attrib
                parameter = child.attrib.get('id')
                if parameter.lower() == 'all':
                    param_list = all_params
                try:
                    width = float(bounds['width'])
                    if len(param_list) > 0:
                        for param_id,value in param_list:
                            val_p = float(value)
                            x = float(width)
                            model.AddResidual(Residuals.PriorInLog('prior_on_%s' % param_id, param_id, scipy.log(val_p),
                                                                   scipy.log(x)))
                            prior_dictionary[param_id] = (val_p, x)
                    else:
                        # TODO: This needs to be fixed
                        param = all_params[all_params.index(parameter)]
                        val_p = param.items()[1]
                        x = width
                        model.AddResidual(Residuals.PriorInLog('prior_on_%s' % param_id, param_id, scipy.log(val_p),
                                                               scipy.log(x)))
                        prior_dictionary[parameter] = (val_p,x)
                except KeyError:
                    lower = float(bounds['lower'])
                    upper = float(bounds['upper'])

                    val_p = sqrt(lower * upper)
                    x = sqrt(upper / val_p)
                    # Just in case something goes wrong but doesn't cause any blatant errors
                    assert (x == sqrt(val_p / lower))
                    model.AddResidual(Residuals.PriorInLog('prior_on_%s' % parameter, parameter, scipy.log(val_p),
                                                           scipy.log(x)))
                    prior_dictionary[parameter] = (val_p,x)
    except KeyError as a:
        logger.warn('Prior Calculation malfunction: ensure prior bounds are numbers, and do not include things like'
                    '"log(x)", "sin(x)", "sqrt()"')
        logger.warn(a)

    return prior_dictionary

# Action routines begin here
# --------------------------


def add_residuals(root, model):
    root = root.find("Parameters")
    prior_dictionary = prior_drier(root, model)




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


def construct_ensemble(params, m, result_dictionary, autocorrelate=False, steps=750, only_pruned=True, prune=1,
                       network_dictionary = None, **kwargs):
    try:
        hess = result_dictionary['hessian']
    except KeyError:
        j = m.jacobian_log_params_sens(log(params))  # use JtJ approximation to Hessian
        hess = dot(transpose(j), j)
    try:
        if kwargs.pop("use_hessian"):
            kwargs['hess'] = hess
    except KeyError:
        pass

    Network.full_speed()
    print 'Beginning ensemble calculation.'
    ens, gs, r = Ensembles.ensemble_log_params(m, params, steps=steps, **kwargs)

    print 'Finished ensemble calculation.'
    if autocorrelate:
        Plotting.figure()
        Plotting.title("Autocorrelation")
        ac = Ensembles.autocorrelation(gs)
        Plotting.plot(ac)
    if prune == 0:
        prune = int(math.ceil(steps/float(250)))
    pruned_ens = ens[::int(prune)]
    if only_pruned:
        if autocorrelate:
            return pruned_ens, ac
        else:
            return pruned_ens
    else:
        return ens, gs, r, pruned_ens

        # print pruned_ens[:, 1]  # finds all values related to r3, goes by index


def plot_histograms(pruned_ens, keyed_list, id='all', bins=10, log=True):
    if id == 'all':
        id = keyed_list.keys()
    if isinstance(id, list):
        for param in id:
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
        param_index = keyed_list.index_by_key(id)
        fig = figure()
        fig.suptitle(id, fontweight='bold')
        param_array = pruned_ens[:, param_index]
        if log:
            hist(scipy.log(param_array), normed=True, bins=bins)
        else:
            hist(param_array, normed=True, bins=bins)


def plot_variables(pruned_ens, net, id=None, time_r=65, start=0, points=200,
                   make_figure=True, color='g', bounds=(.025, .975), net_ensemble=False, graphing_set=None, **kwargs):
    # TODO: make plots show up based on file input.
    times = scipy.linspace(int(start), int(time_r), int(points))
    if net_ensemble:
        bt, mt, st = Ensembles.net_ensemble_trajs(net,
                                                  times,
                                                  pruned_ens)
        vars = []
        for attributes in graphing_set:
            vars.append(attributes["id"])
        Plotting.figure()
        Plotting.plot_ensemble_trajs(bt, mt, st,
                                     vars=vars)
        Plotting.show()
        traj_set = (bt,mt,st)
    else:

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

        print "Plotting trajectory set for %s" % id
        for quantile in quantiles:
            plot(times, quantile.get_var_traj(id), color)

    return traj_set


def calculate_traj(network_dictionary, net, lower_bound = 0, upper_bound = 100):
    net = network_dictionary[net]
    traj = Dynamics.integrate(net, [int(lower_bound), int(upper_bound)])
    return traj, lower_bound, upper_bound


def cost_lm(params, m, optimize=True, plot=False, initial_cost = False, order = None, **kwargs):
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
    optimization_dictionary = {'nelder-mead': Optimization.fmin_log_params,
                               'levenburg-marquardt least square': Optimization.leastsq_log_params,
                               'levenburg-marquardt': Optimization.fmin_lm_log_params}
    return_dictionary = {}
    try:
        plot_after = kwargs.pop('plot_after')
    except KeyError:
        plot_after = False
    if initial_cost:
        initial_cost = m.cost(params)
        print 'Initial Cost:', initial_cost
        try:
            if kwargs.pop('plot_before'):
                initial_plot = Plotting.figure()

                Plotting.plot_model_results(m)
                Plotting.title('Before Optimization')
                return_dictionary["initial"] = initial_plot
        except KeyError:
            pass
    if optimize:
        order.reverse()
        new_params = params.copy()
        for opt in order:
            routine_dict_n = routine_dict_drier(opt.attrib)
            try:
                opt_type = routine_dict_n.pop('type').lower()
            except KeyError:
                opt_type = 'levenburg-marquardt'

            new_params = optimization_dictionary[opt_type](m, new_params, **routine_dict_n)
        optimized_cost = m.cost(new_params)
        params = new_params
        print 'Optimized cost:', optimized_cost
        # print 'Optimized parameters:', params

    if plot or plot_after:
        if not optimize:
            optimized_cost = m.cost(params)
            print 'Optimized cost:', optimized_cost
            print 'Optimized parameters:', params
        Plotting.figure()
        Plotting.plot_model_results(m)
        Plotting.title('After Optimization')
    return_dictionary["params"] = params
    return return_dictionary


def time_extract(experiment_r):
    """
    Goes through the experiment provided and finds the times
    """
    time_array = []
    for experiment in experiment_r.values():

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


def save_to_temp(obj, file_name, xml_file, node, hash_id, routine="temp_file"):
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
    filename = '..\\temp\\' + routine + "_" + str(hash_id) + ".bp"
    Utility.save(obj, filename)
    node.set('path', filename)
    xml_file.write(file_name)


def check_to_save(routine_function):
    """
    Decorator for save/load functionality of action routines
    """

    # TODO: Identify changes to variables and automatically neglect to load from file
    # TODO: Allow multiple files to be chosen from, maybe given names and dynamically accessed from command line?

    def wrapper(root=None, routine=None, hash_node=None, **kwargs):
        """
        Handles the common "path", "save", and "use_path" attributes of most action routines
        
        :param root: (XMLTree) The root of the XML tree
        :param routine: (String) Name of the action routine to be done
        :param kwargs: Extra arguments that will be passed through the decorator function mostly unused
        :return: the result of the associated routine function is returned.  If an object is loaded it is simply
                 passed down to the associated routine function.
        """
        try:
            pass_root = root.find(routine)
            routine_dict = pass_root.attrib
        except AttributeError:
            try:
                current_root = kwargs['current_root']
                pass_root = current_root  # Prevents issue later on where pass_root doesn't work
                routine_dict = pass_root.attrib
            except AttributeError:
                routine_dict = root.attrib
                pass_root = root
        routine_dict = routine_dict_drier(routine_dict)
        use_file = False
        save_file = True
        hash_id = hash_routine(pass_root)
        parent = kwargs['parent']
        hash_root = parent.find('hash')
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
        if hash_root is not None:
            routine_path = hash_root.find(routine.lower())
        else:
            routine_path = None

        if routine_path is not None:
            attributes = routine_path.attrib
            attributes = attributes.copy()
            saved_file_path = attributes.pop('path')
            if os.path.isfile(saved_file_path):
                # A file path is specified
                hash_from_file = os.path.splitext(os.path.basename(saved_file_path))[0].split('_')[1]
                if int(hash_from_file) == hash_id:
                    use_file = True
                    save_file = False
                else:
                    use_file = False
                    save_file = True
                if "override" in sys.argv:
                    # Allows the user to override the cache to avoid saving or loading
                    print "Overriding"
                    use_file = False
                    save_file = False
                elif "replace" in sys.argv:
                    print "replacing"
                    use_file = False
                    save_file = True
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
            action_node = hash_node.find(routine.lower())
            if action_node is None:
                print "Making new node"
                action_node = ET.SubElement(hash_node, routine.lower())
            save_to_temp(routine_object, hash_id=hash_id, xml_file=xml_file,
                         file_name=file_name, node=action_node, routine=routine)
        return routine_object

    return wrapper


@check_to_save
def optimization(routine_dict, model, params, current_root, **kwargs):
    """
    Action routine for optimization.  Checks if the decorator has loaded an object, and if not, runs the regular
    optimization procedure.
    :return: the optimized parameters are returned and passed back to the decorator
    """
    try:
        loaded_object = kwargs['loaded_object']
        optimized_params = loaded_object['params']
        try:
            initial_plot = loaded_object['initial']
            routine_dict.pop("initial_cost")
        except KeyError:
            pass
        return cost_lm(optimized_params, model, optimize=False, **routine_dict)
    except KeyError:
        new_params = params.copy()
        for child in current_root:
            # TODO: Add support for multiple top-level optimizations
            order = recursive_list_builder(child, [child])
            return cost_lm(new_params, model, order=order, **routine_dict)

       # return cost_lm(params, model, **routine_dict)


@check_to_save
def ensemble(routine_dict, result_dictionary, model, network_dictionary, **kwargs):
    """
    Action routine for ensemble construction.  Checks if an ensemble has been loaded, and if not, checks if the 
    optimized parameters have been calculated.
    """
    try:
        # Todo: Find a way to load the autocorrelation graph
        pruned_ensemble = kwargs['loaded_object']
        try:

            Plotting.figure()
            Plotting.title("Autocorrelation")
            Plotting.plot(pruned_ensemble['ac'])
        except IndexError:
            pass
        return pruned_ensemble
    except KeyError:
        pass
    try:
        optimized_params = result_dictionary['optimization']['params']
    except KeyError:
        # We can just pass the unoptimized parameters to the construction function in this case.
        # Eventually should just call its dependent.
        logger.warn("Optimized parameters not found for ensemble construction, using unoptimized parameters instead")
        optimized_params = kwargs['params']
    try:
        autocorrelate = routine_dict.pop('autocorrelate')
        if autocorrelate:
            pruned_ens, ac = construct_ensemble(optimized_params, model, result_dictionary, autocorrelate=autocorrelate,
                                            network_dictionary=network_dictionary, **routine_dict)
            return {'pe': pruned_ens, 'ac': ac}
        else:
            raise KeyError
    except KeyError:
        pruned_ens = construct_ensemble(optimized_params, model, result_dictionary, **routine_dict)
    return {'pe':pruned_ens}


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
        pruned_ens= asarray(pruned_ens['pe'])
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
def ensemble_traj(current_root, routine_dict, result_dictionary, time_r, net, network_dictionary, **kwargs):
    print "Graphing trajectories"
    try:
        traj_dict = kwargs['loaded_object']
        object_loaded = True
    except KeyError:
        object_loaded = False
        traj_dict = {}
    color_array = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    try:
        pruned_ens = result_dictionary['ensemble']['pe']

    except KeyError:
        # We haven't constructed an ensemble yet.
        # TODO: Same as histogram, need ensemble to be called before this in some way
        pass
    try:
        time_r = routine_dict['time']
    except KeyError:
        logger.debug("No time specified, using default")
        # User did not specify a time, so we use the estimated one from the experiment
    # Need to dig into the root to get the variable information
    try:
        points = routine_dict["points"]
    except KeyError:
        points = 100
    try:
        net_to_pass = network_dictionary[routine_dict["net"]]
    except KeyError:
        # No net specified, use default
        net_to_pass = net

    for child_r in current_root:
        for child in child_r:
            graphing_set = []
            index = 0
            try:
                attributes = routine_dict_drier(child.attrib)
                net_ensemble = attributes['net_ensemble']
            except KeyError:
                net_ensemble = False
            for variable in child:
                # We don't run the save function so the attributes don't get dried in the decorator
                # We do it here instead
                attributes = variable.attrib
                attributes = routine_dict_drier(attributes)
                if 'color' not in attributes.keys():
                    attributes['color'] = color_array[index]
                    index += 1
                # We have 8 colors to choose from
                if index > 7:
                    index = 0
                graphing_set.append(attributes)
            figure_bool = True
            if net_ensemble:
                traj_set = plot_variables(pruned_ens, net_to_pass, time_r=time_r, make_figure=figure_bool,
                                          net_ensemble=net_ensemble, graphing_set=graphing_set)
            else:
                for species_attributes in graphing_set:
                    species_id = species_attributes['id']
                    if not object_loaded:
                        traj_set = plot_variables(pruned_ens, net_to_pass, time_r=time_r, make_figure=figure_bool,
                                                  graphing_set=graphing_set, points=points, **species_attributes)
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


@check_to_save
def trajectory_integration(result_dictionary, current_root, network_dictionary, **kwargs):
    traj_list = {}
    for var in current_root.iter("traj"):
        network,lower_bound, upper_bound = calculate_traj(network_dictionary, **var.attrib)
        traj_list[var.attrib['net']] = (network, lower_bound, upper_bound)
    current_column = 0
    for root in current_root.iter("Graph"):
        Plotting.figure(figsize=(15, 4))
        columns = len(root)
        for traj_node in root:
            rows = len(traj_node)
            current_column += 1
            current_row = 0
            for sub in traj_node:
                current_row += 1
                attributes = sub.attrib
                Plotting.subplot(rows, columns, current_column*current_row)
                try:
                    network_id = traj_node.attrib['net']
                except KeyError:
                    network_id = kwargs['net'].id
                traj_and_bounds = traj_list[network_id]
                traj = traj_and_bounds[0]
                if current_row==1:
                    Plotting.title(network_id)

                lower_bound = traj_and_bounds[1]
                upper_bound = traj_and_bounds[2]
                id_list = []
                for var in sub.iter("var"):
                    id_list.append(var.attrib['id'])
                attributes = routine_dict_drier(attributes)
                Plotting.plot_trajectory(traj,id_list, **attributes)
                # Axis solving
                all_vars = traj.dynamicVarKeys
                for thing in traj.assignedVarKeys:
                    all_vars.append(thing)  # The dynamic keys + assigned keys are all the var keys necessary
                max_list = []
                min_list = []
                # Here we find the max and min on the y-axis to constrain the graph
                for var in id_list:
                    index = all_vars.index(var)
                    value_list = []
                    for value in traj.values:
                        value_list.append(value[index])
                    max_list.append(round(max(value_list),1))
                    min_list.append(scipy.floor(min(value_list)))
                padding = max(max_list) * .1  # Pad the max y-value so things look nice
                Plotting.axis([int(lower_bound), int(upper_bound), min(min_list), max(max_list)+padding])


@check_to_save
def hessian(routine_dict, result_dictionary, model, network_dictionary, **kwargs):
    # Todo: Make this load saved object
    try:
        optimized_params = result_dictionary['optimization']['params']
    except KeyError:
        optimized_params = kwargs['params']
    try:
        hess = kwargs['loaded_object']
    except KeyError:
        try:
            eps = routine_dict.pop('eps')
        except KeyError:
            eps = 0.0001
        try:
            if routine_dict.pop('log'):
                hess = model.hessian_log_params(optimized_params, eps, **routine_dict)
        except KeyError:
            hess = model.hessian(optimized_params, eps, **routine_dict)
    evals, evects = Utility.eig(hess)
    Plotting.figure()
    Plotting.title("Eigenvalues")
    Plotting.plot_eigvals(evals)
    Plotting.figure()
    Plotting.title("Eigenvectors")
    Plotting.plot_eigvect(evects[:,0], optimized_params.keys())
    return hess



@check_to_save
def create_Network(current_root, routine_dict, sbml_reference, network_dictionary, network_func_dictionary, **kwargs):
    try:
        network = kwargs['loaded_object']
        return network
    except KeyError:
        pass
    network_name = routine_dict['id']
    if 'from_file' in routine_dict:
        if routine_dict['from_file']:
            network = IO.from_SBML_file(sbml_reference, network_name)
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
                                  'ensembletrajectories': ensemble_traj, 'trajectory': trajectory_integration,
                                  'hessian':hessian}
    network_func_dictionary = {'set_constant': set_constant, 'add_species': add_species,
                               'add_assignment': add_assignment, 'set_optimizable': set_optimizable,
                               'set_initial': set_initial, 'start_fixed': start_fixed,
                               'add_event': add_event, 'add_parameter': add_parameter,
                               'add_rate_rule': add_rate_rule, 'add_reaction': add_reaction}
    # TODO: Clean this mess up
    # These actions should generally happen for any input #
    # They have few enough options that there is little customization to be done and they are essential to setting
    # up the model.


    ###
    hash_node = root.find('hash')
    if hash_node is None:
        hash_node = ET.Element('hash')
        root.append(hash_node)

    # TODO: Build a tree of dependents so that we don't load from file when a dependent changes.
    # TODO: Don't include "independent" operations when calculating hash so that changing them has no effect
    action_root = root.find("Actions")
    network_root = root.find("Networks")
    model_root = root.find("Model")
    net = IO.from_SBML_file(sbml_reference)  # Set base network from SBML reference
    if network_root is not None:
        for network in network_root:
            net_name = network.attrib['id']
            new_network = create_Network(root=network_root, routine=net_name, xml_file=xml_file,
                                                          file_name=file_name, current_root=network,
                                                          sbml_reference=sbml_reference,
                                                          network_dictionary=network_dictionary,
                                                          network_func_dictionary=network_func_dictionary,
                                                          hash_node=hash_node, parent=root)
            network_dictionary[net_name] = new_network

            try:
                attributes = routine_dict_drier(network.attrib)
                if attributes["from_file"]:
                    net = new_network # override default network
            except KeyError:
                pass
    else:
        network_dictionary[net.id] = net
    time_r = time_extract(experiment)

    # TODO: Extend to allow multiple models?
    if model_root is not None:
        experiments = []
        networks = []
        for attribute in model_root:
            if attribute.tag == 'experiment':
                expt_name = attribute.attrib['id']
                expt = experiment[expt_name]
                experiments.append(expt)
            if attribute.tag == 'network':
                net_name = attribute.attrib['id']
                networks.append(network_dictionary[net_name])
        model = Model(experiments, networks)
    else:
        # If undefined, we default to making a model out of all available experiments and models
        model = Model(experiment.values(), network_dictionary.values())

    try:
        fit_root = root.find("Parameters").find("Fit")
        set_ic(fit_root, net)
        params = key_parameters(fit_root)
    except AttributeError:
        # No parameters were specified, so we just take them all
        params = model.get_params()
    # This function is standalone and just tacks a residual onto anything that needs it


    add_residuals(root, model)


    for child in action_root:
        #try:
        r_name = child.tag.lower()
        result_dictionary[r_name] = action_function_dictionary[r_name](root=action_root, routine=child.tag,
                                                                    model=model, params=params, current_root=child,
                                                                    xml_file=xml_file, file_name=file_name, time_r=time_r,
                                                                    net=net, result_dictionary=result_dictionary,
                                                                    parent=root, hash_node=hash_node,
                                                                    network_dictionary = network_dictionary)
        #except KeyError as e:
            #logger.warn("No function associated with action tag %s" % child.tag)
            #logger.warn(e)
    show()


            #
        # plot_histograms(pruned, ['r3', 'tao'], params)
        # show()
