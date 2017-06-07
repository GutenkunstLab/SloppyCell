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
from scipy import *
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


def construct_ensemble(params, m, num_steps=750, only_pruned=True):
    j = m.jacobian_log_params_sens(log(params))  # use JtJ approximation to Hessian
    jtj = dot(transpose(j), j)
    print 'Beginning ensemble calculation.'
    ens, gs, r = Ensembles.ensemble_log_params(m, asarray(params), jtj, steps=num_steps)
    print 'Finished ensemble calculation.'

    pruned_ens = asarray(ens[::int(math.ceil(num_steps/float(250)))])
    if only_pruned:
        return pruned_ens
    else:
        return ens, pruned_ens

        # print pruned_ens[:, 1]  # finds all values related to r3, goes by index


def plot_histograms(pruned_ens, params, keyed_list, bins=10):
    if isinstance(params, list):
        for param in params:
            if param in keyed_list.keys():
                param_index = keyed_list.index_by_key(param)
                fig = figure()
                fig.suptitle(param, fontweight='bold')
                param_array = pruned_ens[:, param_index]
                hist(log(param_array), normed=True, bins=bins)
    else:
        param_index = keyed_list.index_by_key(params)
        fig = figure()
        fig.suptitle(params , fontweight='bold')
        param_array = pruned_ens[:, param_index]
        hist(log(param_array), normed=True, bins=bins)


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


def cost_lm(params, m, optimize=True, plotter=False):

    if optimize:
        initial_cost = m.cost(params)
        print 'Initial Cost:', initial_cost
        params = Optimization.fmin_lm_log_params(m, params, maxiter=20, disp=False)
        optimized_cost = m.cost(params)
        print 'Optimized cost:', optimized_cost
        print 'Optimized parameters:', params
    if plotter:
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


def find_by_tag(root, tag, attr=True, return_node = False):
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
    print "Saving " + routine + " to file"
    filename = '..\\temp\\' + routine + "_" + str(int(math.ceil(time.time()))) + ".bp"
    Utility.save(obj, filename)
    node.set('path', filename)
    xml_file.write(file_name)


def check_to_save(routine_function):
    """
    Decorator for save functionality of action routines
    """
    def wrapper(root=None, routine=None, **kwargs):
        """
        Handles the common "path", "save", and "use_path" attributes of most action routines
        
        :param root: (XMLTree) The root of the XML tree
        :param routine: (String) Name of the action routine to be done
        :param kwargs: Extra arguments that will be passed through the decorator function mostly unused
        :return: the result of the associated routine function is returned
        """
        routine_dict = root.find(routine).attrib
        use_file = False
        save_file = False

        try:
            if routine_dict['use_path'].lower() == 'true':
                use_file = True
        except KeyError as X:
            logger.info("Path use not specified")
            logger.debug(X)
        try:
            if routine_dict['save'].lower() == 'true':
                save_file = True
        except KeyError as X:
            logger.info("Save attribute not specified, will only save if no path currently exists.")
            logger.debug(X)
        if 'path' in routine_dict:
            if os.path.isfile(routine_dict['path']):
                # A file path is specified
                if use_file:
                    # User wants to load the object
                    print 'Successfully loaded %s objects from file' % routine
                    loaded_object = Utility.load(routine_dict['path'])
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
                            os.unlink(routine_dict['path'])
                            print "Replacing old %s file with new one" % routine
                        except KeyError:
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
    try:
        plotter = routine_dict['plot']
        if plotter == 'True':
            plotter = True
        else:
            plotter = False
    except KeyError:
        plotter = False
    try:
        optimized_params = kwargs['loaded_object']
        return cost_lm(optimized_params, model, optimize=False, plotter=plotter)
    except KeyError:
        return cost_lm(params, model, plotter=plotter)


@check_to_save
def ensemble(routine_dict, result_dictionary, model, **kwargs):
    try:
        pruned_ensemble = kwargs['loaded_object']
        return pruned_ensemble
    except KeyError:
        pass
    try:
        optimized_params = result_dictionary['optimization']
    except KeyError:
        logger.warn("Optimized parameters not found for ensemble construction, using unoptimized parameters instead")
        optimized_params = kwargs['params']
    try:
        num_steps = int(routine_dict['size'])
        return construct_ensemble(optimized_params, model, num_steps)
    except KeyError:
        return construct_ensemble(optimized_params, model)


def make_happen(root, experiment, xml_file=None, file_name=None):
    result_dictionary = dict()
    function_dictionary = {'optimization': optimization, 'ensemble': ensemble}
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
                                                                    model=model, params=params,
                                                                    xml_file=xml_file, file_name=file_name,
                                                                    result_dictionary=result_dictionary)
        except KeyError:
            logger.warn("No function associated with action tag %s" % child.tag)
    show()

        #
        # plot_histograms(pruned, ['r3', 'tao'], params)
        # show()
