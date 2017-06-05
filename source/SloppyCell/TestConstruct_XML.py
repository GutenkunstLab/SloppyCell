"""
Created 6/2/2017

Author @Keeyan

Runs appropriate functions to test input data and model (XML edition
"""

import re
from pylab import *
from scipy import *
from SloppyCell.ReactionNetworks import *
from xml.etree import ElementTree as ET
import logging

logger = logging.getLogger('TestConstruct_XML')
logging.basicConfig()

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
    pruned_ens = asarray(ens[::len(params) * 5])
    if only_pruned:
        return pruned_ens
    else:
        return ens, pruned_ens

        # print pruned_ens[:, 1]  # finds all values related to r3, goes by index


def plot_histograms(pruned_ens, params, keyed_list):
    if isinstance(params, list):
        for param in params:
            if param in keyed_list.keys():
                param_index = keyed_list.index_by_key(param)
                fig = figure()
                fig.suptitle(param, fontweight='bold')
                param_array = pruned_ens[:, param_index]
                hist(log(param_array), normed=True)
    else:
        param_index = keyed_list.index_by_key(params)
        fig = figure()
        fig.suptitle(params , fontweight='bold')
        param_array = pruned_ens[:, param_index]
        hist(log(param_array), normed=True)


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


def cost_lm(params, m, optimize=True, plotter=True):
    print 'Initial Cost:', m.cost(params)
    if optimize:
        params = Optimization.fmin_lm_log_params(m, params, maxiter=20, disp=False)
        print 'Optimized cost:', m.cost(params)
        print 'Optimized parameters:', params
    if plotter:
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


def make_happen(root, experiment):
    sbml_reference = root.find("References").find("SBML").attrib["path"]
    net = IO.from_SBML_file(sbml_reference)  # Set model from SBML reference
    time = time_extract(experiment)
    model = Model([experiment], [net])  # Use experiment to create model
    fit_root = root.find("Parameters").find("Fit")
    set_ic(fit_root, net)
    add_residuals(root, model)
    params = key_parameters(fit_root)
    optimized_params = cost_lm(params, model)
    pruned = construct_ensemble(optimized_params, model)
    plot_histograms(pruned, ['r3', 'tao'], params)
    show()
