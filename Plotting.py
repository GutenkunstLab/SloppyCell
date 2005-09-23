import scipy
from pylab import *
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

eigvals_cW = ColorWheel(colors = ('b', 'g', 'c', 'm', 'k'), 
                        lines = None)
def plot_eigvals(vals, label=None, offset = 0):
    posVals = abs(scipy.compress(vals > 0, vals))
    posRange = scipy.compress(vals > 0, range(len(vals)))
    negVals = abs(scipy.compress(vals < 0, vals))
    negRange = scipy.compress(vals < 0, range(len(vals)))

    sym = eigvals_cW.next()
    line = semilogy(posRange+offset, posVals, sym, label = label)
    if len(negVals) > 0:
        semilogy(negRange+offset, negVals, 'r'+sym[1:])

    a = axis()
    axis([-0.05*len(vals) + offset, 1.05*(len(vals) - 1) + offset, a[2], a[3]])

    return line

PlotEigenvalueSpectrum = plot_eigvals

def plot_eigvect(vect, labels=None, num_label = 5):
    """
    Plot a given eigenvector.

    If a list of labels is passed in, the largest (in magnitude) num_label bars
     will be labeled on the plot.
    """
    # The 0.4 centers the bars on their numbers, accounting for the default
    #  bar width of 0.8
    bar(scipy.arange(len(vect)) - 0.4, vect/scipy.linalg.norm(vect))
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
