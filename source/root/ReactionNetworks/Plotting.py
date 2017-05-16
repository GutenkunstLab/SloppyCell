import sets

import scipy

import SloppyCell.Plotting
from SloppyCell.Plotting import *
import Network_mod

import logging
logger = logging.getLogger('ReactionNetworks.Plotting')

def PlotEigenvectors(eigVects, net = None, title = None):
    nEv = 3
    nOv = len(eigVects[:,0])
    for jj in range(nEv):
        subplot(nEv, 1, jj+1)
        if jj == 0 and title is not None:
            title(title)

        bar(range(nOv), eigVects[:,jj]/scipy.linalg.norm(eigVects[:,jj]))
        axis([-1, nOv] + axis()[2:])

        if net is not None:
            mags = zip(abs(eigVects[:,jj]), range(nOv), eigVects[:,jj])
            mags.sort()
            mags.reverse()
            for mag, index, val in mags[:5]:
                name = net.optimizableVars[index].name
                if name is None:
                    name = net.optimizableVars[index].id
                text(index, val + scipy.sign(val)*0.05, 
                           name,
                           horizontalalignment='center',
                           verticalalignment='center')

        a = list(axis())
        a[0:2] = [-.03*nOv, nOv*1.03]
        a[2] -= 0.1
        a[3] += 0.1
        axis(a)

def plotStateSpaceTrajectoriesForVariables(traj, id1, id2, thresholds = None):
    xx = traj.getVariableTrajectory(id1)
    yy = traj.getVariableTrajectory(id2)
    plot(xx, yy)
    plot([xx[0]], [yy[0]], 'or')
    xlabel(id1)
    ylabel(id2)
    if thresholds is not None:
        a = axis()
        vlines([thresholds[0]], a[2], a[3])
        hlines([thresholds[1]], a[0], a[1])

def plotTrajectoriesForVariables(traj, ids = None, showLegend = True):
    if ids is None:
        ids = traj.net.variables.keys()

    cW = ColorWheel()

    lines = []
    legend = []
    for id in ids:
        line = plot(traj.timepoints, traj.getVariableTrajectory(id), 
                          cW.next()[::2])
        lines.append(line)
        legend.append(id)

    if showLegend:
        legend(lines, legend)


def PlotTrajectoriesForExperiments(model, experiments, params = None, with_data=True,
                                   plotPts = 100, overlap = .1, skip = 1, showLegend=True):

    # First find the maximum time in our data
    maxTime = 0
    chemsNeededByCalc = {}
    for exptName in experiments:
        dataByCalc = model.exptColl[exptName].GetData()
        for calc in dataByCalc:
            chemsNeededByCalc.setdefault(calc, [])
            for chem in dataByCalc[calc].keys():
                chemsNeededByCalc[calc].append(chem)
                thisMaxTime = max(dataByCalc[calc][chem].keys()) 
                if thisMaxTime > maxTime:
                    maxTime = thisMaxTime

    lines = []
    legend = []
    times = scipy.linspace(0, maxTime*(1 + overlap), plotPts)
    varsByCalc = {}
    for calc in chemsNeededByCalc:
        varsByCalc[calc] = {}
        for chem in chemsNeededByCalc[calc]:
            varsByCalc[calc][chem] = times

    model.GetCalculationCollection().Calculate(varsByCalc, params)
    calcVals = model.GetCalculationCollection().GetResults(varsByCalc)
    cW = ColorWheel()
    for exptName in experiments:
        expt = exptColl[exptName]
        dataByCalc = expt.GetData()
        for calc in dataByCalc:
            for chem in dataByCalc[calc]:
                color, sym, dash = cW.next()

                if with_data:
                    for time, (data, error) in dataByCalc[calc][chem].items()[::skip]:
                        errorbar(time, data, yerr=error, color=color,
                                 mfc = color, marker=sym,
                                 ecolor=color, capsize=6)

                predicted = scipy.array(calcVals[calc][chem].items())
                order = scipy.argsort(predicted[:,0])
                predicted = scipy.take(predicted, order)
                predicted[:,1] = predicted[:,1] *\
                        model.GetScaleFactors()[exptName][chem]
                lines.append(plot(predicted[:,0], predicted[:,1],
                                  color=color, linestyle=dash,
                                  linewidth = 3))
                legend.append(chem + ' in ' + str(calc))# + ' for ' + str(exptName))


    if showLegend:
        legend(lines, legend, loc=4)

def PlotDataForExperiments(model, experiments, skip = 1):
    exptColl = model.GetExperimentCollection()

    cW = ColorWheel()
    for exptName in experiments:
        expt = exptColl[exptName]
        dataByCalc = expt.GetData()
        for calc in dataByCalc:
            for chem in dataByCalc[calc]:
                color, sym, dash = cW.next()
                d = scipy.zeros((len(dataByCalc[calc][chem].values()[::skip]),
                                            3), scipy.Float)
                for ii, (time, (data, error))\
                        in enumerate(dataByCalc[calc][chem].items()[::skip]):
                    d[ii] = [time, data, error]

                errorbar(d[:,0], d[:,1], yerr=d[:,2], color=color,
                         mfc=color, marker=sym,
                         ecolor=color, capsize=6)

def plot_model_data(model, expts = None, style = 'errorbars',
                    show_legend = True, loc = 'upper left'):
    """
    Plot the data in the given experiments for the given model.
    
    Note: You may need to run a Plotting.show() to display the plot.
    
    Inputs:
      model: Model whose experiments to plot
      expts: List of experiment IDs to plot
      style: Style of plot. Currently supported options are:
          'errorbars': Plots points and bars for each data point
          'lines': Plots a continuous line for the data
      show_legend: Boolean that control whether or not to show the legend
      loc: Location of the legend. See help(Plotting.legend) for options.
    """
    return plot_model_results(model, expts, style, show_legend, loc, 
                              plot_trajectories = False)

def plot_model_results(model, expts = None, style='errorbars',
                       show_legend = True, loc = 'upper left',
                       plot_data = True, plot_trajectories = True,
                       data_to_plot=None):
    """
    Plot the fits to the given experiments for the last cost evalution of the
    model.
    
    Note: You may need to run a Plotting.show() to display the plot.
    
    Inputs:
      model: Model whose results to plot
      expts: List of experiment IDs to plot, if None is specified, all 
             experiments are plotted
      style: Style of plot. Currently supported options are:
              'errorbars': Plots points and bars for each data point
              'lines': Plots a continuous line for the data
      show_legend: Boolean that control whether or not to show the legend
      loc: Location of the legend. See help(Plotting.legend) for options.
      plot_data: Boolean that controls whether the data is plotted
      plot_trajectories: Boolean that controls whether the trajectories are
                         plotted
      data_to_plot: If None, all data variables will be plotted. Otherwise,
                       pass a list of id's and only variables in  that list
                       will be plotted.
    """
    exptColl = model.get_expts()
    calcColl = model.get_calcs()

    lines, labels = [], []
    cW = ColorWheel()
    
    if expts is None:
        expts = exptColl.keys()

    for exptId in expts:
        expt = exptColl[exptId]
        dataByCalc = expt.GetData()
        for ds in expt.scaled_extrema_data:
            dataByCalc.setdefault(ds['calcKey'], {})
            dataByCalc[ds['calcKey']].setdefault(ds['var'], {})
        # We sort the calculation names for easier comparison across plots
        sortedCalcIds = dataByCalc.keys()
        sortedCalcIds.sort()
        for calcId in sortedCalcIds:
            # Pull the trajectory from that calculation, defaulting to None
            #  if it doesn't exist.
            net = calcColl.get(calcId)
            traj = getattr(net, 'trajectory', None)
            for dataId, dataDict in dataByCalc[calcId].items():
                # Skip this variable if it's not amongst our data_to_plot
                #  list.
                if (data_to_plot is not None) and (dataId not in data_to_plot):
                    continue
                color, sym, dash = cW.next()

                if plot_trajectories:
                    if traj is None:
                        print 'No trajectory in calculation %s!' % calcId
                        print 'The cost must be evaluated before the results',
                        print 'can be plotted.'
                        return
                    scaleFactor = model.GetScaleFactors()[exptId][dataId]
                    result = scaleFactor*traj.getVariableTrajectory(dataId)
                    l = plot(traj.timepoints, result, color=color, 
                             linestyle=dash,linewidth=3)

                    # We superimpose a dotted black line to distinguish
                    #  theory from data in this case
                    if style is 'lines':
                        plot(traj.timepoints, result, 'k--', linewidth=3,
                                   zorder = 10)
                

                if plot_data and dataDict:
                    # Pull the data out of the dictionary and into an array
                    d = scipy.array([[t, v, e] for (t, (v, e))
                                     in dataDict.items()])
                    if style is 'errorbars':
                        l = errorbar(d[:,0], d[:,1], yerr=d[:,2], color=color,
                                     mfc = color, linestyle='', marker=sym,
                                     ecolor='k', capsize=6)[0]
                    elif style is 'lines':
                        # Make sure we order the data before plotting
                        order = scipy.argsort(d[:,0], 0)
                        d = scipy.take(d, order, 0)
                        l = plot(d[:,0], d[:,1], color=color,
                                 linestyle=dash)
                        lines.append(l)

                # Plot the extra data points.
                if plot_data:
                    for res in model.residuals:
                        if isinstance(res, Residuals.ScaledExtremum)\
                           and res.exptKey == exptId and res.calcKey == calcId\
                           and res.var == dataId:
                            t = res.last_time_result
                            val = res.yMeas
                            sigma = res.ySigma
                            errorbar([t], [val], [sigma], color=color,
                                     linestyle='', marker=sym,
                                     ecolor=color, capsize=6, mfc='w', 
                                     mec=color, mew=2, ms=10)

                # Let's print the pretty name for our variable if we can.
                lines.append(l)
                name = net.get_component_name(dataId)
                labels.append('%s in %s for %s' % (name, calcId, exptId))

    if show_legend:
        legend(lines, labels, loc=loc)

    return lines, labels

def plot_ensemble_results(model, ensemble, expts = None, 
                          style='errorbars',
                          show_legend = True, loc = 'upper left',
                          plot_data = True, plot_trajectories = True):
    """
    Plot the fits to the given experiments over an ensemble. 

    Note that this recalculates the cost for every member of the ensemble, so
     it may be very slow. Filtering correlated members from the ensemble is
     strongly recommended.
    
    Inputs:
      model: Model whose results to plot
      ensemble: Parameter ensemble
      expts: List of experiment IDs to plot, if None is specified, all 
             experiments are plotted
      style: Style of plot. Currently supported options are:
              'errorbars': Plots points and bars for each data point
              'lines': Plots a continuous line for the data
      show_legend: Boolean that control whether or not to show the legend
      loc: Location of the legend. See help(Plotting.legend) for options.
      plot_data: Boolean that controls whether the data is plotted
      plot_trajectories: Boolean that controls whether the trajectories are
                         plotted
    """
    exptColl = model.get_expts()
    nets = model.get_calcs()

    if expts is None:
        expts = exptColl.keys()

    lines, labels = [], []
    cW = ColorWheel()

    Network_mod.Network.pretty_plotting()
    model.cost(ensemble[0])
    timepoints = {}
    for netId, net in nets.items():
        traj = getattr(net, 'trajectory', None)
        if traj is not None:
            net.times_to_add = scipy.linspace(traj.timepoints[0], 
                                              traj.timepoints[-1], 1000)

    Network_mod.Network.full_speed()

    results = {}
    for params in ensemble:
        model.cost(params)
        for exptId in expts:
            expt = exptColl[exptId]
            results.setdefault(exptId, {})
            dataByCalc = expt.GetData()
            for netId in dataByCalc.keys():
                results[exptId].setdefault(netId, {})
                # Pull the trajectory from that calculation, defaulting to None
                #  if it doesn't exist.
                net = nets.get(netId)
                traj = net.trajectory
                for dataId in dataByCalc[netId].keys():
                    results[exptId][netId].setdefault(dataId, [])

                    scaleFactor = model.GetScaleFactors()[exptId][dataId]
                    result = scaleFactor*traj.get_var_traj(dataId)
                    results[exptId][netId][dataId].append(result)

    for exptId in expts:
        expt = exptColl[exptId]
        dataByCalc = expt.GetData()
        # We sort the calculation names for easier comparison across plots
        sortedCalcIds = dataByCalc.keys()
        sortedCalcIds.sort()
        for netId in sortedCalcIds:
            for dataId, dataDict in dataByCalc[netId].items():
                color, sym, dash = cW.next()

                if plot_data:
                    # Pull the data out of the dictionary and into an array
                    d = scipy.array([[t, v, e] for (t, (v, e))
                                     in dataDict.items()])
                    if style is 'errorbars':
                        l = errorbar(d[:,0], d[:,1], yerr=d[:,2], fmt='o', 
                                     color=color, markerfacecolor=color,
                                     marker=sym, ecolor='k', capsize=6)[0]
                    elif style is 'lines':
                        # Make sure we order the data before plotting
                        order = scipy.argsort(d[:,0], 0)
                        d = scipy.take(d, order, 0)
                        l = plot(d[:,0], d[:,1], color=color,
                                 linestyle=dash)
		    lines.append(l)

                if plot_trajectories:
                    times = model.get_calcs().get(netId).trajectory.get_times()
                    mean_vals = scipy.mean(results[exptId][netId][dataId], 0)
                    std_vals = scipy.std(results[exptId][netId][dataId], 0)

                    lower_vals = mean_vals - std_vals
                    upper_vals = mean_vals + std_vals

                    # Plot the polygon
                    xpts = scipy.concatenate((times, times[::-1]))
                    ypts = scipy.concatenate((lower_vals, upper_vals[::-1]))
                    fill(xpts, ypts, fc=color, alpha=0.4)

                # Let's print the pretty name for our variable if we can.
                name = net.get_component_name(dataId)
                labels.append('%s in %s for %s' % (name, netId, exptId))

    for netId, net in nets.items():
        del net.times_to_add

    if show_legend:
        legend(lines, labels, loc=loc)

    for net in nets.values():
        net.times_to_add = None

    return lines, labels


def plot_trajectory(traj, vars = None,
                    show_legend = True, loc = 'upper left',
                    logx = False, logy = False):
    if vars is None:
        vars = traj.dynamicVarKeys

    plot_funcs_dict = {(False, False): plot,
                       (True, False): semilogx,
                       (False, True): semilogy,
                       (True, True): loglog}
    plot_func = plot_funcs_dict[(logx, logy)]

    cW = ColorWheel(symbols=None)
    lines, labels = [], []
    for id in vars:
        color, sym, dash = cW.next()
        label = str(id)
        line = plot_func(traj.timepoints, traj.getVariableTrajectory(id),
                         color = color, mfc = color, marker=sym, linestyle=dash,
                         linewidth=3, label=label)
        lines.append(line)
        labels.append(label)

    if show_legend:
        legend(loc=loc)

    return (lines, labels)

def plot_deriv_trajectory(traj, vars = None,
                    show_legend = True, loc = 'upper left',
                    logx = False, logy = False):
    if vars is None:
        vars = traj.dynamicVarKeys

    plot_funcs_dict = {(False, False): plot,
                       (True, False): semilogx,
                       (False, True): semilogy,
                       (True, True): loglog}
    plot_func = plot_funcs_dict[(logx, logy)]

    cW = ColorWheel(symbols=None)
    lines, labels = [], []
    for id in vars:
        try:
            color, sym, dash = cW.next()
            label = str(id)
            line = plot_func(traj.timepoints, traj.get_var_traj((id, 'time')),
                             color = color, mfc = color, marker=sym, linestyle=dash,
                             linewidth=3, label=label)
            lines.append(line)
            labels.append(label)
        except ValueError:
            logger.warn('Derivative of variable %s not present in trajectory. Make sure '\
                        'return_derivs is set to True in integrate function.'  % (id))
            

    if show_legend:
        legend(loc=loc)

    return (lines, labels)

def plot_ensemble_trajs(best_traj=None, mean_traj=None, 
                        std_traj=None, std_devs=1.0,
                        vars=None, 
                        show_legend = True, loc = 'upper left'):
    """
    Plot the results of a Network ensemble.

    Inputs:
     best_traj -- Best-fit trajectory
     mean_traj -- Mean trajectory
     std_traj -- Trajectory of standard deviations
     std_devs -- Number of standard deviations to draw bounds at
     vars -- List of variable ids to plot the bounds for

     show_legend -- Boolean to show a legend or not
     loc -- Location code for legend
    """
    cW = ColorWheel()
    lines = []
    labels = []
    for var in vars:
        color, sym, dash = cW.next()
        if best_traj is not None:
            l = plot(best_traj.timepoints, best_traj.getVariableTrajectory(var),
                     linestyle='-', color=color, linewidth=2)
        if mean_traj is not None:
            times = mean_traj.timepoints
            mean_vals = mean_traj.getVariableTrajectory(var)
            plot(times, mean_vals, linestyle = '--', linewidth=2,
                 color=color)
        if mean_traj is not None and std_traj is not None and std_devs > 0:
            # Calculate the bounds
            std_vals = std_traj.getVariableTrajectory(var)
            lower_vals = mean_vals - std_devs*std_vals
            upper_vals = mean_vals + std_devs*std_vals

            # Plot the polygon
            xpts = scipy.concatenate((times, times[::-1]))
            ypts = scipy.concatenate((lower_vals, upper_vals[::-1]))
            fill(xpts, ypts, facecolor=color, alpha=0.4)

        lines.append(l)
        labels.append(var)

    if show_legend:
        legend(lines, labels, loc=loc)
