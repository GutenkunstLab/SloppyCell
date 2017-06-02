"""
Created 5/30/2017

Author @Keeyan

Runs appropriate functions to test input data and model
"""


import re
from pylab import *
from scipy import *
from SloppyCell.ReactionNetworks import *
import logging
logger = logging.getLogger('TestConstruct')
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


def prior_drier(dictop):
    """
    Finds the appropriate values to put into a PriorInLog function to constrain the function.  
    Simple algebra tells us that 'val_p' is equal to the square root of the upper bound multiplied
    by the lower bound, and that 'x' is equal to the square root of the upper bound divided by 
    'val_p', or the square root of the 'val_p' divided by the lower bound.  
    The log of 'val_p' and 'x' become the parameters necessary for the PriorInLog function to bound
    the parameters correctly (the p value and sigma p value, respectively)
    
    :param dictop: Full dictionary output from ScellParser, contains the prior dictionary
    :return: A dictionary mapping each parameter to its respective (val_p,x) pair.  The (val_p,x) pair is 
      contained in a dictionary for convenience 
    """
    # TODO: Maybe don't use regex?  Allow for log and other functions to be used in prior definition
    priors = dictop['priors']
    # This nifty regex (sorry) works in the following way:
    # '[^,]' causes it to not match leading commas (the ones that separate values and aren't located inside parentheses)
    # '.*?:' matches to any sequence of characters followed by a colon, the '.*?' expression
    # means "any amount of characters"
    # '\(.*?\)' the backslashes escape special characters.  The '\(.*?\) means "any amount of characters in between
    # two parentheses.  The question mark ensures its the shortest string it can find, so that we don't
    # accidentally count other parentheses later on as part of '.*', which matches any character.
    # Without the question mark it would match the full string.
    entries = re.findall('[^,].*?:\(.*?\)', priors)
    # An entry has the format 'parameter_name:(lower_bound,upper_bound)'
    prior_dictionary = {}
    val_p = 0
    x = 0
    for entry in entries:
        parameter = entry.split(':')[0]
        try:
            constraint_list = []
            # The regex expression basically chops off the leading and trailing parentheses by matching for everything
            # except a leading and trailing parentheses, using '[^(]' and '[^)]'.  Regex was probably not needed.
            # Chopping off the first and last character by index seemed too barbaric.
            boundaries = re.findall('[^(].*[^)]', entry.split(':')[1])[0].split(',')

            for bound in boundaries:
                constraint_list.append(float(bound))
            lower = min(constraint_list)  # just in case the bounds are in the wrong order for whatever reason
            upper = max(constraint_list)
            val_p = sqrt(lower*upper)
            x = sqrt(upper/val_p)
            assert(x == sqrt(val_p/lower))  # Just in case something goes wrong but doesn't cause any syntactical errors
        except Exception as a:
            logger.warn('Prior Calculation malfunction: ensure prior bounds are numbers, and do not include things like'
                        '"log(x)", "sin(x)", "sqrt()"')
            logger.warn(a)
        # Don't necessarily have to make a dictionary, but it makes things easier to read.
        prior_dictionary[parameter] = {'val': val_p, 'x': x}
    return prior_dictionary


def add_residuals(full_dict, model):
    prior_dictionary = prior_drier(full_dict)
    for param in prior_dictionary.keys():
        values = prior_dictionary[param]
        res = Residuals.PriorInLog(param+'_prior', param, log(values['val']), log(values['x']))
        model.AddResidual(res)


def set_ic(param_list, net):
    for param in param_list:
        net.set_var_ic(chop(param, '_0'), param)


def find_vars(param_list, characters):
    prior_list = []
    for param in param_list:
        if param.lower().endswith(characters):
            prior_list.append(chop(param, characters))
    return prior_list


def construct_ensemble(params, m, num_steps=750, only_pruned=True):
    j = m.jacobian_log_params_sens(log(params))  # use JtJ approximation to Hessian
    jtj = dot(transpose(j), j)
    print 'Beginning ensemble calculation.'
    ens, gs, r = Ensembles.ensemble_log_params(m, asarray(params), jtj, steps=num_steps)
    print 'Finished ensemble calculation.'
    pruned_ens = asarray(ens[::len(params)*5])
    if only_pruned:
        return pruned_ens
    else:
        return ens, pruned_ens

    # print pruned_ens[:, 1]  # finds all values related to r3, goes by index



def plot_histograms(pruned_ens, param_index):
    figure()
    hist(log(pruned_ens[:, param_index]), normed=True)
    pass


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


def time_extract(untagged_dict):
    """
    Goes through the experiment provided and figures out how long 
    """
    time_array = []
    for key_t in untagged_dict['experiment'].data.keys():
        for species_key in untagged_dict['experiment'].data[key_t].keys():
            times = untagged_dict['experiment'].data[key_t][species_key].keys()
            time_array.append(max(times))
    return max(time_array) + len(untagged_dict['parameters_to_fit'])


def make_happen(untagged_dict):
    net = IO.from_SBML_file(untagged_dict['SBML_Reference'])  # Set model from SBML reference
    time = time_extract(untagged_dict)
    model = Model([untagged_dict['experiment']], [net])  # Use experiment to create model
    fit_list = untagged_dict['parameters_to_fit']
    set_ic(find_vars(fit_list.keys(), '*'), net)
    add_residuals(untagged_dict, model)
    param_list = []
    # Place parameters into a list of tuples to be place into a KeyedList
    for key_p in untagged_dict['parameters_to_fit'].keys():
        chopped_key = chop(key_p, ['*', '^'])
        param_list.append((chopped_key, untagged_dict['parameters_to_fit'].get(key_p)))
    params = KeyedList(param_list)

    optimized_params = cost_lm(params, model)
    construct_ensemble(optimized_params, model, net, time=time)
