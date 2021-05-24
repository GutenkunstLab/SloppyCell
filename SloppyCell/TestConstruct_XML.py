"""
Created 6/2/2017

Author @Keeyan

Runs appropriate functions to test input data and model (XML edition
"""
from __future__ import print_function
from __future__ import absolute_import

import logging
from xml.etree import ElementTree as ET
import os
import scipy
import time
from SloppyCell.ReactionNetworks import *
from pylab import *
import csv

from . import Utility

logger = logging.getLogger('TestConstruct_XML')
logging.basicConfig()
this_module = sys.modules[__name__]
# Globally defined variables.  Should eventually be done away with.
step_factor = 300
saved_already = []


def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def findoutlier(a, sensitivity=3.5):
    third_quartile = scipy.percentile(a, 75)
    first_quartile = scipy.percentile(a, 25)
    iqr = third_quartile-first_quartile
    outlier_list = []
    for number in a:
        if number > (third_quartile+iqr*sensitivity) or number < (first_quartile-iqr*sensitivity):
            outlier_list.append(number)
    return outlier_list


def find_factors(x):
    n = math.sqrt(x)
    n = math.floor(n)
    while x % n != 0:
        n -= 1
    m = x/n
    return m, n


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
        return recursive_list_builder(child, order)


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
                        # Here we take care of any conditionals that are used for network creation
                        if key_r == "condition":
                            a = ['!=', '==', '<', '>', '*', '/']
                            matches = next((x for x in a if x in parameter), False)
                            value = float(parameter.split(matches)[1])
                            b = {'!=': lambda x: x != value, '==': lambda x: x == value, '<': lambda x: x < value,
                                 '>': lambda x: x > value, '*': lambda x: x*value, '/': lambda x: x/value}
                            parameter = b[matches]

        routine_dict_copy[key_r] = parameter
    return routine_dict_copy

# Here we define all of the modification functions for the networks
# ------------------------------------------------------------------------


def add_compartment(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.add_compartment(**attributes)
    return network


def set_dynamic(network, action_root, network_dictionary):
    for child in action_root:
        if child.tag == 'traj':
            attributes = child.attrib
            traj, up_b, low_b = calculate_traj(network_dictionary, **attributes)
            vals = traj.last_dynamic_var_values()
            network.set_dyn_var_ics(vals)
    return network


def set_typical(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        if attributes['id'] == 'all':
            variables = network.variables.keys()
            for variable in variables:
                network.set_var_typical_val(variable, attributes['value'])
        else:
            network.set_var_typical_val(**attributes)

    return network


def set_constant(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.set_var_constant(**attributes)
    return network


def add_species(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.add_species(**attributes)
    return network


def add_assignment(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
        var_id = attributes.pop('id')
        network.add_assignment_rule(var_id, **attributes)
    return network


def set_optimizable(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        network.set_var_optimizable(**attributes)
    return network


def set_initial(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        try:
            condition = attributes.pop('condition')
        except KeyError:
            condition = lambda x: True
        if attributes['id'] == 'all':
            variables = network.variables.keys()
            for variable in variables:
                alpha = condition(network.get_var_ic(variable))
                if isinstance(alpha,bool):
                    if alpha:
                        network.set_var_ic(variable, attributes['value'])
                else:
                    network.set_var_ic(variable, alpha)
        elif attributes['id']=='traj' and 'point' in attributes.keys():
            point = int(attributes['point'])
            traj=Dynamics.integrate(network,[0,point])
            network.set_var_ics(traj.get_var_vals(point))
        else:
            if 'value' in attributes:
                if condition(attributes['value']):
                    network.set_initial_var_value(**attributes)
            else:
                value = network.get_var_ic(attributes['id'])
                alpha = condition(value)
                if isinstance(alpha,bool):
                    network.set_var_ic(value=value, **attributes)
                else:
                    network.set_var_ic(value=alpha,**attributes)
    return network


def start_fixed(network, action_root, network_dictionary):
    attributes = action_root.attrib
    attributes = routine_dict_drier(attributes)
    try:
        min_r = int(attributes['min'])
    except KeyError:
        min_r = 0
    try:
        max_r = int(attributes['max'])
    except KeyError:
        max_r = 1000
    Dynamics.integrate(network, [min_r, max_r])
    fp = Dynamics.dyn_var_fixed_point(network)
    network.set_dyn_var_ics(fp)
    return network


def add_event(network, action_root, network_dictionary):
    # TODO: This could be easier to use
    for child in action_root:
        attributes = child.attrib
        assignment_dict = {}
        for assignment in child:
            attribs = assignment.attrib
            assignment_dict[attribs['id']] = attribs['func']
        network.add_event(event_assignments=assignment_dict, **attributes)
    return network


def add_parameter(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = routine_dict_drier(attributes)
        try:
            listed_network = attributes.pop('net')
            listed_network = network_dictionary[listed_network]
            variables = listed_network.parameters.keys()
            for id in variables:
                network.add_parameter(id, listed_network.get_var_ic(id), **attributes)
        except KeyError:
            network.add_parameter(**attributes)
    return network


def add_rate_rule(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
        var_id = attributes.pop('id')
        network.add_rate_rule(var_id, **attributes)
    return network


def add_reaction(network, action_root, network_dictionary):
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
            pass_list.append(stoich)
        try:
            reaction_name = attributes.pop('reaction')
            try:
                kinetic_law = getattr(globals()['Reactions'], reaction_name)
                network.addReaction(kinetic_law, var_id, *pass_list, **attributes)
                return network
            except Exception as exep:
                print(exep)
                logger.warn("Please check that the reaction name is typed correctly")
        except KeyError:
            # No reaction specified
            kinetic_law = attributes.pop('kinetic_law')
            network.addReaction(var_id, kinetic_law, *pass_list, **attributes)
            return network
        for key_r in attributes:
            pass_list.append(attributes[key_r])
        network.addReaction(var_id, *pass_list)
    return network


def remove_component(network, action_root, network_dictionary):
    for child in action_root:
        attributes = child.attrib
        attributes = attributes.copy()
        attributes = routine_dict_drier(attributes)
        network.remove_component(**attributes)
        return network

# -----------------------------------------------------------------


def prior_drier(root, model):
    """
    Finds the appropriate values to put into a PriorInLog function to constrain the function.  
    Simple algebra tells us that 'val_p' is equal to the square root of the upper bound multiplied
    by the lower bound, and that 'x' is equal to the square root of the upper bound divided by 
    'val_p', or the square root of the 'val_p' divided by the lower bound.  
    The log of 'val_p' and 'x' become the parameters necessary for the PriorInLog function to bound
    the parameters correctly (the p value and sigma p value, respectively)

    :param root: Full dictionary output from ScellParser, contains the prior dictionary
    :param model: The model that we want to add priors to
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
                        for param_id, value in param_list:
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
                    # assert (x == sqrt(val_p / lower))
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


def plot_traj(traj_dict, network_dictionary, id_list, attributes, id=None,
              log=False, lower_bound=0, upper_bound=100, points=1, file_name=None, network_id=None, vary=True):
    if vary:
        if log:
            space = scipy.logspace(int(lower_bound), int(upper_bound), int(points))
        else:
            space = scipy.linspace(int(lower_bound), int(upper_bound), int(points))
    else:
        space = [1]
    network = network_dictionary[traj_dict['net']]
    if len(id_list)<1:
        id_list = None

    for point in space:
        if vary:
            network.set_var_ic(id,point)
        traj, lower_bound, upper_bound = calculate_traj(network_dictionary, **traj_dict)
        try:
            semilog = attributes["semilog"]
            agg = semilogger(traj, id_list, semilog)
            a = Plotting.semilogx(traj.get_times(), agg)
            semilogged = True
        except KeyError:
            semilogged = False
            agg=None
            a = Plotting.plot_trajectory(traj, id_list, **attributes)
    if file_name is not None:
        dir_name = os.path.dirname(file_name)
        parent_name = os.path.basename(dir_name)
        hash_id = 0
        hash_id = id_list[0]+"_"+str(upper_bound)+"_"+str(len(id_list))
        traj.to_file(dir_name+"/saved_files/Traj_"+network_id+"_"+hash_id+".csv")
    return traj, lower_bound, upper_bound, semilogged, a, agg



def semilogger(traj,ids, semilog):
    id_copy = list(ids)
    agg=traj.get_var_traj(id_copy.pop(0))
    for id in ids:
        agg=agg+traj.get_var_traj(id)
    agg = 1-agg
    return agg



def add_residuals(root, model):
    root = root.find("Parameters")
    if root is not None:
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
                       network_dictionary = None, file_name = None, **kwargs):
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
    print('Beginning ensemble calculation.')
    print(params)
    ens, gs, r = Ensembles.ensemble_log_params(m, params, steps=steps, **kwargs)

    print('Finished ensemble calculation.')
    if autocorrelate:
        print(autocorrelate)
        Plotting.figure()
        Plotting.title("Autocorrelation")
        ac = Ensembles.autocorrelation(gs)
        Plotting.plot(ac)
    else:
        ac = None
    if prune == 0:
        prune = int(math.ceil(steps/float(250)))
    pruned_ens = ens[::int(prune)]
    if file_name is not None:
        dir_name = os.path.dirname(file_name)
        csv_friendly_ens = asarray(pruned_ens)
        with open(dir_name+"/saved_files/Ensemble_ens_"+str(steps)+".csv", 'wb') as myfile:
            wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
            print(params.keys())
            wr.writerow(params.keys())
            for row in csv_friendly_ens:
                wr.writerow(row)
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
                print("Plotting histogram for %s" % param)
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


def plot_variables(pruned_ens=None, net=None, id=None, time_r=65, start=0, points=200,
                   make_figure=True, color='g', bounds=(.025, .975), net_ensemble=False, graphing_set=None, **kwargs):
    # TODO: make plots show up based on file input.
    times = scipy.linspace(int(start), int(time_r), int(points))
    if net_ensemble:
        bt, mt, st = Ensembles.net_ensemble_trajs(net,
                                                  times,
                                                  pruned_ens)
        # We could write the bt, mt, and st to file
        # bt.to_file("henlo.csv")
        vars = []
        for attributes in graphing_set:
            vars.append(attributes["id"])
        Plotting.figure()
        Plotting.plot_ensemble_trajs(bt, mt, st,
                                     vars=vars)
        Plotting.title("Ensemble Trajectory")
        traj_set = (bt, mt, st)
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

        print("Plotting trajectory set for %s" % id)
        for quantile in quantiles:
            plot(times, quantile.get_var_traj(id), color)
    return traj_set


def calculate_traj(network_dictionary, net, lower_bound = 0, upper_bound = 100, reset=False):
    net = network_dictionary[net]
    if reset:
        net.resetDynamicVariables()
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
        print(('Initial Cost:', initial_cost))
        try:
            if kwargs.pop('plot_before'):
                initial_plot = Plotting.figure()

                f=Plotting.plot_model_results(m)
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
        print(('Optimized cost:', optimized_cost))
        # print 'Optimized parameters:', params

    if plot or plot_after:
        if not optimize:
            optimized_cost = m.cost(params)
            print(('Optimized cost:', optimized_cost))
            print(('Optimized parameters:', params))
        Plotting.figure()
        f=Plotting.plot_model_results(m)
        for thing in f[0]:
            thing.set_alpha(.7)

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
    root = xml_file.getroot()
    try:
        output = root.find("saved_files").attrib["path"]
    except KeyError:
        output = None
    if output is not None:
        dir_name = output
    else:
        dir_name = os.path.dirname(file_name)
    parent_name = os.path.basename(dir_name)
    print("Saving " + routine + " to file")
    print("Output folder: " + dir_name)
    model_name = root.attrib['name']
    save_folder = "/saved_files/" + routine +"-"+model_name+ "_" +str(hash_id)+ ".bp"

    filename = dir_name+save_folder

    relative_path=dir_name+save_folder
    if not os.path.isdir(dir_name+"/saved_images"):
        os.mkdir(dir_name+"/saved_images")
    for i in plt.get_fignums():
        if i not in saved_already:
            saved_already.append(i)
            Plotting.figure(i)
            save_image_folder = '/saved_images/'+routine+"-"+model_name+'-figure%d.png' % i
            Plotting.savefig(dir_name+save_image_folder)
    try:
        Utility.save(obj, filename)
    except IOError:
        os.mkdir(dir_name+"/saved_files/")
        Utility.save(obj, filename)
    node.set('path', relative_path)
    root = xml_file.getroot()
    indent(root)
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
        hash_root = parent.find('saved_files')
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
                hash_from_file = os.path.splitext(os.path.basename(saved_file_path))[0].split('_')[-1]
                if int(hash_from_file) == hash_id:

                    use_file = True
                    save_file = False
                else:
                    use_file = False
                    save_file = True
                if "override" in sys.argv:
                    # Allows the user to override the cache to avoid saving or loading
                    print("Overriding")
                    use_file = False
                    save_file = False
                elif "replace" in sys.argv:
                    print("replacing")
                    use_file = False
                    save_file = True
                if use_file:
                    # User wants to load the object
                    print('Successfully loaded %s objects from file' % routine)
                    print(saved_file_path)
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
                            print("Replacing old %s file with new one" % routine)
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
                print("Making new node")
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
    print(type(params))
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
def ensemble(routine_dict, result_dictionary, model, network_dictionary, file_name, **kwargs):
    """
    Action routine for ensemble construction.  Checks if an ensemble has been loaded, and if not, checks if the 
    optimized parameters have been calculated.
    """
    try:
        # Todo: Find a way to load the autocorrelation graph
        pruned_ensemble = kwargs['loaded_object']
        try:
            ac=pruned_ensemble['ac']
            Plotting.figure()
            Plotting.title("Autocorrelation")
            Plotting.plot(ac)
        except KeyError:
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
                                                network_dictionary=network_dictionary, file_name=file_name,
                                                **routine_dict)
            return {'pe': pruned_ens, 'ac': ac}
        else:
            raise KeyError
    except KeyError:
        pruned_ens = construct_ensemble(optimized_params, model, result_dictionary, file_name=file_name, **routine_dict)
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
    print("Graphing trajectories")
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
                        plot_variables(traj_set=traj_set, make_figure=figure_bool, time_r=time_r, points=points,
                                       **species_attributes)

                    if figure_bool:
                        figure_bool = False
    return traj_dict


@check_to_save
def trajectory_integration(result_dictionary, current_root, network_dictionary, file_name=None, **kwargs):
    # TODO: Break up into smaller functions
    try:
        traj_list = kwargs['loaded_object']
    except KeyError:
        traj_list = {}

    color_array = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    color_pick = 0
    current_column = 0
    traj_count=0
    for root in current_root.iter("Graph"):
        columns = 0
        sub_row_array = []
        for traj_node in root.iter("traj"):
            traj_count+=1
            sub_rows = 0
            column_list = []
            for thing in traj_node:
                column_list.append(thing.tag)
                if thing.tag == "subplot":
                    sub_rows += 1
            if "subplot" in column_list:
                columns = columns + 1
            sub_row_array.append(sub_rows)
        sub_rows = min(sub_row_array)
        if sub_rows == 0:
            sub_rows = 1
        current_column = 0
        if columns == 0:
            columns = 1
        if columns == 1:
            sub_rows, sub_columns = find_factors(sub_rows)
            current_fig = Plotting.figure(figsize=(sub_columns * 5, sub_rows * 4))
        elif sub_rows == 1:
            sub_columns, sub_rows = find_factors(columns)
            current_fig = Plotting.figure(figsize=(sub_columns * 4, sub_rows * 5))
        else:
            current_fig = Plotting.figure(figsize=(columns*5,sub_rows*4))
        max_list = []
        min_list = []
        placeholder_columns = columns
        for traj_node in root.iter("traj"):
            columns=placeholder_columns
            rows = len(traj_node)
            if columns == 1:
                rows, columns = find_factors(rows)
            elif rows == 1:
                if sub_rows>columns:
                    rows, columns = find_factors(placeholder_columns)
                else:
                    columns, rows = find_factors(placeholder_columns)
            column_list = []
            for thing in traj_node:
                column_list.append(thing.tag)
            if "subplot" in column_list:
                current_column = current_column + 1
            current_row = 1
            run_once = True
            indicator=False
            for sub in traj_node:
                id_list = []
                if sub.tag.lower() == "subplot":

                    color_pick = 0
                    max_list = []
                    min_list = []

                    attributes = sub.attrib
                    if current_row>rows*columns:
                        if current_column==1:
                            current_row=current_row-(rows*columns)+1
                        elif current_row > rows and current_row * current_column > rows * columns:
                            current_row = current_row - (rows)
                    elif current_row*current_column>rows*columns:
                        current_row=current_row-(columns)+1
                    Plotting.subplot(rows, columns, current_column*current_row)
                    current_row += columns
                    for var in sub.iter("var"):
                        id_list.append(var.attrib['id'])
                else:
                    attributes = root.attrib
                    indicator=True
                    for var in traj_node.iter("var"):
                        id_list.append(var.attrib['id'])
                if run_once:
                    if indicator:
                        run_once = False
                    try:
                        network_id = traj_node.attrib['net']
                    except KeyError:
                        network_id = kwargs['net'].id
                    varied = False
                    attributes = routine_dict_drier(attributes)
                    try:
                        show_outlier = attributes.pop("show_outlier")
                    except KeyError:
                        show_outlier = False
                    try:
                        lower_margin = attributes.pop("lower_bound")
                    except KeyError:
                        lower_margin = None
                    try:
                        upper_margin = attributes.pop("upper_bound")
                    except KeyError:
                        upper_margin = None
                    traj_dict = routine_dict_drier(traj_node.attrib)
                    try:
                        style=traj_dict.pop('style')
                        color=traj_dict.pop('color')
                    except KeyError:
                        style=None
                        color=None
                    for vary in sub.iter("vary"):
                        varied = True
                        traj, lower_bound, upper_bound, semilogged, a, agg = plot_traj(traj_dict, network_dictionary,
                                                                                       id_list, attributes,
                                                                                       network_id=network_id,file_name=file_name, **vary.attrib)
                    if not varied:
                        traj, lower_bound, upper_bound, semilogged, a, agg = plot_traj(traj_dict, network_dictionary,
                                                                                       id_list, attributes,
                                                                                       network_id=network_id, file_name=file_name, vary=False)

                    if current_row == 1+columns:
                        Plotting.title(network_id)

                    if "upper_bound" in traj_dict:
                        upper_bound = int(traj_dict['upper_bound'])



                    # Plotting.plot(traj.get_times(),
                    # traj.get_var_traj(id_list[0]), **attributes)

                    if isinstance(a[0],list):
                        lines = a[0]
                        if len(lines) < 2:
                            linex = lines[0][0]
                            linex.set_color(color_array[color_pick])
                        for line in lines:
                            line=line[0]
                            if style is not None:
                                line.set_linestyle(style)
                            if color is not None:
                                line.set_color(color)

                        legend(loc="upper left")
                        color_pick += 1
                    else:
                        pass
                    #TODO: Make sure lines are separate colors regardless of scenario
                    # Axis solving
                    all_vars = traj.dynamicVarKeys
                    if len(id_list) < 1:
                        id_list = traj.dynamicVarKeys
                    for thing in traj.assignedVarKeys:
                        all_vars.append(thing)  # The dynamic keys + assigned keys are all the var keys necessary

                    # Here we find the max and min on the y-axis to constrain the graph
                    if semilogged:
                        max_list=agg
                        min_list=agg
                    else:
                        for var in id_list:
                            index = all_vars.index(var)
                            value_list = []
                            for value in traj.values:
                                value_list.append(value[index])
                            max_list.append(round(max(value_list), 3))
                            if min(value_list)>0:
                                min_list.append(scipy.floor(min(value_list)))
                            else:
                                min_list.append(round(min(value_list),3))


                    if not show_outlier:
                        while len(findoutlier(max_list)) > 0:
                            for outlier in findoutlier(max_list):
                                max_list.remove(outlier)
                    padding = max(max_list) * .1  # Pad the max y-value so things look nice
                    if upper_margin is None:
                        upper_margin=max(max_list)+padding
                    if lower_margin is None:
                        lower_margin=min(min_list)


                    Plotting.axis([float(lower_bound), float(upper_bound), lower_margin, upper_margin])
        current_fig.subplots_adjust(hspace=.3)
    return traj_list


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
    network = None
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
        network = network_func_dictionary[r_name](network, child, network_dictionary)
    if "compile" in routine_dict:
        network.compile()


    return network

def scale_factors(scale_root,expts):
    print("Scaling")
    print("\n")
    experiments = {}
    for expt in expts:
        experiments[expt.name] = expt
    for expt in scale_root:
        current_expt = experiments[expt.attrib['id']]
        scale_dict = {}
        shared_dict = {}
        for variable in expt:
            attrib = variable.attrib
            if('value' in attrib):
                scale_dict[variable.attrib['id']]=float(variable.attrib['value'])
            if('shared' in attrib):
                try:
                    shared_dict[attrib['shared']].append(attrib['id'])
                except Exception as e:
                    shared_dict[attrib['shared']] = []
                    shared_dict[attrib['shared']].append(attrib['id'])
        current_expt.set_shared_sf(shared_dict.values())
        current_expt.set_fixed_sf(scale_dict)
def make_happen(root, experiment, xml_file=None, file_name=None, sbml_reference = None, output_location=None):
    result_dictionary = dict()
    network_dictionary = dict()
    file_name = os.path.abspath(file_name)
    action_function_dictionary = {'optimization': optimization, 'ensemble': ensemble, 'histogram': histogram_r,
                                  'ensembletrajectories': ensemble_traj, 'trajectory': trajectory_integration,
                                  'hessian': hessian}
    network_func_dictionary = {'set_constant': set_constant, 'add_species': add_species,
                               'add_assignment': add_assignment, 'set_optimizable': set_optimizable,
                               'set_initial': set_initial, 'start_fixed': start_fixed,
                               'add_event': add_event, 'add_parameter': add_parameter,
                               'add_rate_rule': add_rate_rule, 'add_reaction': add_reaction,
                               'set_typical': set_typical, "remove": remove_component, 'set_dynamic': set_dynamic,
                               'add_compartment': add_compartment}
    # TODO: Clean this mess up
    # These actions should generally happen for any input #
    # They have few enough options that there is little customization to be done and they are essential to setting
    # up the model.
    ###
    hash_node = root.find('saved_files')
    if hash_node is None:
        hash_node = ET.Element('saved_files')
        root.append(hash_node)
    if output_location is not "unspecified" and output_location is not None:
        hash_node.set("path",output_location)
    # TODO: Build a tree of dependents so that we don't load from file when a dependent changes.
    # TODO: Don't include "independent" operations when calculating hash so that changing them has no effect
    action_root = root.find("Actions")
    network_root = root.find("Networks")
    model_root = root.find("Model")
    end_session = False
    if network_root is None:
        logger.critical("A network section is required for use of XML format input.  Use the <Network> tag and set its ID to 'base' to default to the sbml file")
        end_session = True
    if model_root is None:
        end_session =True
        logger.critical("A model section is required for use of XML format.  Use the <Model> tag and include an experiment and at least one network")
    if not end_session:

        try:
            net = IO.from_SBML_file(sbml_reference)  # Set base network from SBML reference
        except ValueError:
            net = None
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
        if experiment is not None:
            time_r = time_extract(experiment)
        else:
            time_r = 0
        # TODO: Extend to allow multiple models?
        try:
            get_from_model=False
            fit_root = root.find("Parameters").find("Fit")
            params = key_parameters(fit_root)
            for network in network_dictionary:
                network = network_dictionary[network]
                set_ic(fit_root, network)
                for param in network.parameters.keys():
                    if param not in params.keys():
                        network.set_var_constant(param, False)

        except AttributeError:
            get_from_model = True
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
            # If undefined, we default to making a model out of all available experiments and networks
            if experiment is not None:
                model = Model(experiment.values(), network_dictionary.values())
            else:
                # We can only make a model with an experiment defined
                model = None

        if get_from_model is True:
            # No parameters were specified, so we just take them all
            if model is not None:
                params = model.get_params()
            else:
                params = None
        # This function is standalone and just tacks a residual onto anything that needs it
        if model is not None:
            add_residuals(root, model)
        scale_root = root.find("Scale_Factors")
        if(scale_root is not None):
            scale_factors(scale_root,experiments)
        if action_root is not None:

            for child in action_root:
                # try:
                r_name = child.tag.lower()
                action_function = action_function_dictionary[r_name]
                result_dictionary[r_name] = action_function(root=action_root, routine=child.tag,
                                                            model=model, params=params, current_root=child,
                                                            xml_file=xml_file, file_name=file_name, time_r=time_r,
                                                            net=net, result_dictionary=result_dictionary,
                                                            parent=root, hash_node=hash_node,
                                                            network_dictionary = network_dictionary)
                # except KeyError as e:
                #     logger.warn("No function associated with action tag %s" % child.tag)
                #     logger.warn(e)
        print("All routines complete")

    time.sleep(2)
