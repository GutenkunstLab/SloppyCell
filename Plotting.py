import pylab, scipy

def ColorWheel():
    """
    ColorWheel()

    Returns a generator that cycles through a selection of colors, symbols, and 
    line styles for matlibplot.matlab.plot.
    """
    colors = ['b', 'g', 'r', 'c', 'm', 'k']
    symbols = ['o', 's', '^', 'v', '<', ">", 'x', 'D', 'h', 'p']
    lines = ['-', '--', '-.', ':']
    while 1:
        for l in lines:
           for s in symbols:
                for c in colors:
                   yield c + s + l

cW = ColorWheel()
def PlotEigenvalueSpectrum(vals, lab = None):
    posVals = abs(scipy.compress(vals > 0, vals))
    posRange = scipy.compress(vals > 0, range(len(vals)))
    negVals = abs(scipy.compress(vals < 0, vals))
    negRange = scipy.compress(vals < 0, range(len(vals)))

    sym = cW.next()
    while sym[0] == 'r':
        sym = cW.next()

    line = pylab.semilogy(posRange, posVals, sym[:2], label = lab)
    if len(negVals) > 0:
        pylab.semilogy(negRange, negVals, 'r'+sym[1])

    a = pylab.axis()
    pylab.axis([-3, len(vals)+2, a[2], a[3]])

    return line
