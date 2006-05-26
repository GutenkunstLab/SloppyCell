import scipy
try:
    from pylab import *
except RuntimeError:
    # When running in parallel we found that this import could raise a 
    # 'RuntimeError: could not open display' rather than an ImportError, so
    # we catch and raise an error we know how to handle
    raise ImportError

import Residuals

rc('lines', linewidth=2)

def ColorWheel(colors = ('b', 'g', 'r', 'c', 'm', 'k'), 
               symbols = ('o', 's', '^', 'v', '<', ">", 'x', 'D', 'h', 'p'),
               lines = ('-', '--', '-.', ':')):
    """
    ColorWheel()

    Returns a generator that cycles through a selection of colors, symbols, and 
    line styles for matlibplot.matlab.plot.
    """
    if not colors:
        colors = ('',)
    if not symbols:
        symbols = ('',)
    if not lines:
        lines = ('',)

    while 1:
        for l in lines:
           for s in symbols:
                for c in colors:
                   yield c + s + l

vals_cW = 0
def reset_vals_cw():
    """
    Reset the ColorWheel used for plotting eigenvalues.
    """
    global vals_cW
    vals_cW = ColorWheel(colors = ('b', 'r', 'g', 'c', 'm', 'y', 'k'), 
                         lines = None)
reset_vals_cw()

def plot_eigvals(vals, label=None, offset=0, indicate_neg=True, join=False):
    posVals = abs(scipy.compress(scipy.real(vals) > 0, vals))
    posRange = scipy.compress(scipy.real(vals) > 0, range(len(vals)))
    negVals = abs(scipy.compress(scipy.real(vals) < 0, vals))
    negRange = scipy.compress(scipy.real(vals) < 0, range(len(vals)))

    sym = vals_cW.next()
    if indicate_neg:
        if len(negVals) > 0:
            semilogy(negRange+offset, negVals, 'r'+sym[1:], zorder=1)
        if sym[0] == 'r':
            sym = vals_cW.next()

    line = semilogy(posRange+offset, posVals, sym, label = label, zorder=0)

    if join:
        plot(scipy.arange(len(vals)) + offset, abs(vals), sym[0]+'--', 
             zorder=-1)

    a = axis()
    axis([-0.05*len(vals) + offset, 1.05*(len(vals) - 1) + offset, a[2], a[3]])

    return line

def plot_singvals(vals, label=None, offset=0, join=False):
    return plot_eigvals(vals, label, offset, indicate_neg=False, join=join)

PlotEigenvalueSpectrum = plot_eigvals

def plot_eigvect(vect, labels=None, bottom = 0, num_label = 5):
    """
    Plot a given eigenvector.

    If a list of labels is passed in, the largest (in magnitude) num_label bars
     will be labeled on the plot.
    """
    # The 0.4 centers the bars on their numbers, accounting for the default
    #  bar width of 0.8
    vect = scipy.real(vect)
    max_index = scipy.argmax(abs(vect))
    if vect[max_index] < 0:
        vect = -vect
    bar(scipy.arange(len(vect)) - 0.4, vect/scipy.linalg.norm(vect), 
        bottom=bottom)
    a = axis()
    a[0:2] = [-.03*len(vect) - 0.4, (len(vect) - 1)*1.03 + 0.4]

    if labels is not None:
        mags = zip(abs(vect), range(len(vect)), vect)
        mags.sort()
        mags.reverse()
        for mag, index, val in mags[:num_label]:
            name = labels[index]
            text(index, val + scipy.sign(val)*0.05, name,
                 horizontalalignment='center', verticalalignment='center')

        a[2] -= 0.1
        a[3] += 0.1

    axis(a)

def plot_priors(model,priorIDs=None,params=None,sameScale=False):
    """
    Plots specified priors and parameter values.
    If no priors are specified, plots them all.
    If no params are provided, uses the model params.
    If sameScale is true, recenters everything so all prior optima are at 0.
    Labeling is awkward and hence avoided here.  I suggest using the
    pylab.text command with parameter names after the plot has been generated.
    """
    if params is None:
        params=model.get_params()
    residuals = model.GetResiduals()
    if priorIDs is None:
        priorIDs = residuals.keys()

    priorVals=[]
    priorErrs=[]
    parVals=[]
    for resID in priorIDs:
        res = residuals.getByKey(resID)
        if isinstance(res, Residuals.PriorInLog):
            priorVals.append(res.logPVal)
            priorErrs.append(res.sigmaLogPVal)
            parVals.append(params.getByKey(res.pKey))

    if sameScale is False:
        errorbar(scipy.arange(len(priorVals)),priorVals,yerr=priorErrs,fmt='bo',ecolor='k',capsize=6)
        errorbar(scipy.arange(len(priorVals)),scipy.log(parVals),fmt='go')
    else:
        errorbar(scipy.arange(len(priorVals)),scipy.zeros(len(priorVals)),yerr=priorErrs,fmt=None,ecolor='k',capsize=6)
        errorbar(scipy.arange(len(priorVals)),scipy.log(parVals)-priorVals,fmt='go')

            
