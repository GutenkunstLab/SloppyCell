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
                            print "Seperating Values"
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


def plot_variables(pruned_ens, net, time=65, start=0, points=100):
    # TODO: make plots show up based on file input.
    times = linspace(start, time, points)
    traj_set = Ensembles.ensemble_trajs(net, times, pruned_ens)
    lower, upper = Ensembles.traj_ensemble_quantiles(traj_set, (0.025, 0.975))

    figure()
    plot(times, lower.get_var_traj('frac_v3'), 'g')
    plot(times, upper.get_var_traj('frac_v3'), 'g')
    plot(times, lower.get_var_traj('frac_v4'), 'b')
    plot(times, upper.get_var_traj('frac_v4'), 'b')

    show()


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
    return max(time_array) + 5  # Arbitrary number added, still unsure how x-axis is calculated for plots


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
        routine_dict = root.find(routine).attrib
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
    # Graph objects could potentially be saved, but it's probably not worth it.
    try:
        pruned_ens = result_dictionary['ensemble']
    except KeyError:
        # We might just not have done ensemble construction yet.
        # TODO: Either organize actions before calling, or have actions call their dependents
        pass
    for child in current_root:
        for variable in child:
            attributes = variable.attrib
            attributes = routine_dict_drier(attributes)

            plot_histograms(pruned_ens, params, **attributes)


def make_happen(root, experiment, xml_file=None, file_name=None):
    result_dictionary = dict()
    function_dictionary = {'optimization': optimization, 'ensemble': ensemble, 'histogram': histogram_r}
    # TODO: Clean this mess up
    # These actions should generally happen for any input #
    # They have few enough options that there is little customization to be done and they are essential to setting
    # up the model.
    sbml_reference = root.find("References").find("SBML").attrib["path"]
    net = IO.from_SBML_file(sbml_reference)  # Set model from SBML reference
    time = time_extract(experiment)
    model = Model([experiment], [net])  # Use experiment to create model
    fit_root = root.find("Parameters").find("Fit")
    set_ic(fit_root, net)
    add_residuals(root, model)
    params = key_parameters(fit_root)
    ###
    action_root = root.find("Actions")
    for child in action_root:
        try:
            r_name = child.tag.lower()
            result_dictionary[r_name] = function_dictionary[r_name](root=action_root, routine=child.tag,
                                                                    model=model, params=params, current_root=child,
                                                                    xml_file=xml_file, file_name=file_name,
                                                                    result_dictionary=result_dictionary)
        except KeyError:
            logger.warn("No function associated with action tag %s" % child.tag)
    show()

        #
        # plot_histograms(pruned, ['r3', 'tao'], params)
        # show()
