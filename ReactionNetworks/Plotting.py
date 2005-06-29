
import scipy
from pylab import *

import SloppyCell.Plotting
PlotEigenvalueSpectrum = SloppyCell.Plotting.PlotEigenvalueSpectrum

def PlotEigenvectors(eigVects, net = None, title = None):
    nEv = 3
    nOv = len(eigVects[:,0])
    for jj in range(nEv):
        pylab.subplot(nEv, 1, jj+1)
        if jj == 0 and title is not None:
            pylab.title(title)

        pylab.bar(range(nOv), eigVects[:,jj]/scipy.linalg.norm(eigVects[:,jj]))
        pylab.axis([-1, nOv] + pylab.axis()[2:])

        if net is not None:
            mags = zip(abs(eigVects[:,jj]), range(nOv), eigVects[:,jj])
            mags.sort()
            mags.reverse()
            for mag, index, val in mags[:5]:
                name = net.optimizableVars[index].name
                if name is None:
                    name = net.optimizableVars[index].id
                pylab.text(index, val + scipy.sign(val)*0.05, 
                           name,
                           horizontalalignment='center',
                           verticalalignment='center')

        a = pylab.axis()
        a[0:2] = [-.03*nOv, nOv*1.03]
        a[2] -= 0.1
        a[3] += 0.1
        pylab.axis(a)

def plotStateSpaceTrajectoriesForVariables(traj, id1, id2, thresholds = None):
    xx = traj.getVariableTrajectory(id1)
    yy = traj.getVariableTrajectory(id2)
    pylab.plot(xx, yy)
    pylab.plot([xx[0]], [yy[0]], 'or')
    pylab.xlabel(id1)
    pylab.ylabel(id2)
    if thresholds is not None:
        a = pylab.axis()
        pylab.vlines([thresholds[0]], a[2], a[3])
        pylab.hlines([thresholds[1]], a[0], a[1])

def plotTrajectoriesForVariables(traj, ids = None, showLegend = True):
    if ids is None:
        ids = traj.net.variables.keys()

    cW = SloppyCell.Plotting.ColorWheel()

    lines = []
    legend = []
    for id in ids:
        line = pylab.plot(traj.timepoints, traj.getVariableTrajectory(id), 
                          cW.next()[::2])
        lines.append(line)
        legend.append(id)

    if showLegend:
        pylab.legend(lines, legend)


def PlotTrajectoriesForExperiments(model, experiments, params = None, with_data=True,
                                   plotPts = 100, overlap = .1, skip = 1, showLegend=True):

    # First find the maximum time in our data
    maxTime = 0
    chemsNeededByCalc = {}
    for exptName in experiments:
        dataByCalc = exptColl[exptName].GetData()
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
    cW = SloppyCell.Plotting.ColorWheel()
    for exptName in experiments:
        expt = exptColl[exptName]
        dataByCalc = expt.GetData()
        for calc in dataByCalc:
            for chem in dataByCalc[calc]:
                fmt = cW.next()
                if with_data:
                    for time, (data, error) in dataByCalc[calc][chem].items()[::skip]:
                        pylab.errorbar(time, data, yerr=error, fmt=fmt,
                                         ecolor=fmt, capsize=6)

                predicted = scipy.array(calcVals[calc][chem].items())
                order = scipy.argsort(predicted[:,0])
                predicted = scipy.take(predicted, order)
                predicted[:,1] = predicted[:,1] *\
                        model.GetScaleFactors()[exptName][chem]
                lines.append(pylab.plot(predicted[:,0], predicted[:,1],
                                          fmt[::2], linewidth = 3))
                legend.append(chem + ' in ' + str(calc))# + ' for ' + str(exptName))


    if showLegend:
        pylab.legend(lines, legend, loc=4)

def PlotDataForExperiments(model, experiments, skip = 1):
    exptColl = model.GetExperimentCollection()

    cW = SloppyCell.Plotting.ColorWheel()
    for exptName in experiments:
        expt = exptColl[exptName]
        dataByCalc = expt.GetData()
        for calc in dataByCalc:
            for chem in dataByCalc[calc]:
                fmt = cW.next()
                d = scipy.zeros((len(dataByCalc[calc][chem].values()[::skip]),
                                            3), scipy.Float)
                for ii, (time, (data, error))\
                        in enumerate(dataByCalc[calc][chem].items()[::skip]):
                    d[ii] = [time, data, error]

                pylab.errorbar(d[:,0], d[:,1], yerr=d[:,2], fmt=fmt[:-1],
                                 ecolor=fmt, capsize=6)

def plot_model_data(model, expts = None, style = 'errorbars',
                    show_legend = True, loc = 'upper left'):
    plot_model_results(model, expts, style, show_legend, loc, 
                       plot_trajectories = False)

def plot_model_results(model, expts = None, style='errorbars',
                       show_legend = True, loc = 'upper left',
                       plot_data = True, plot_trajectories = True):
    exptColl = model.GetExperimentCollection()
    calcColl = model.GetCalculationCollection()

    lines, labels = [], []
    cW = SloppyCell.Plotting.ColorWheel()
    
    if expts is None:
        expts = exptColl.keys()

    for exptId in expts:
        expt = exptColl[exptId]
        dataByCalc = expt.GetData()
        sortedCalcIds = dataByCalc.keys()
        sortedCalcIds.sort()
        for calcId in sortedCalcIds:
            traj = getattr(calcColl[calcId], 'trajectory', None)
            for dataId, dataDict in dataByCalc[calcId].items():
                fmt = cW.next()

                if plot_data:
                    # Pull the data out of the dictionary and into an array
                    d = scipy.array([[t, v, e] for (t, (v, e))
                                     in dataDict.items()])
                    if style is 'errorbars':
                        l = errorbar(d[:,0], d[:,1], yerr=d[:,2], 
                                           fmt=fmt[:-1], ecolor=fmt, capsize=6)
                    elif style is 'lines':
                        order = scipy.argsort(d[:,0], 0)
                        d = scipy.take(d, order, 0)
                        l = plot(d[:,0], d[:,1], fmt[::2])

                if plot_trajectories:
                    if traj is None:
                        print 'No trajectory for calculation %s!' % calcId
                        print 'The cost must be evaluated before the results',
                        print 'can be plotted.'
                        return

                    scaleFactor = model.GetScaleFactors()[exptId][dataId]
                    result = scaleFactor*traj.getVariableTrajectory(dataId)
                    plot(traj.timepoints, result, fmt[::2], linewidth=3)

                    if style is 'lines':
                        plot(traj.timepoints, result, 'k--', linewidth=3,
                                   zorder = 10)
                
                lines.append(l)
                labels.append('%s in %s for %s' % (dataId, calcId, exptId))

    if show_legend:
        legend(lines, labels, loc=loc)
