"""
Alters SciPy 0.3.2 to look like version SciPy 0.5.1 as far as SloppyCell is
concerned.

See also addAssignmentRulesToFunctionBody in ReactionNetworks/Network_mod.py
"""
import scipy

if int(scipy.__version__.split('.')[1]) < 4:
    # Type names have changed
    scipy.float_ = scipy.Float
    scipy.float64 = scipy.Float64
    scipy.int_ = scipy.Int
    scipy.int32 = scipy.Int32

    # Create an empty class to move around the limits module
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52308
    class Bunch(dict):
        def __init__(self,**kw):
            dict.__init__(self,kw)
            self.__dict__.update(kw)

    # Move scipy.limits to scipy.misc.limits
    scipy.misc = Bunch(limits = scipy.limits)

    # The seeding for the random number generator has changed. This
    #  should be compatible and useful
    scipy.stats.old_seed = scipy.stats.seed
    def new_seed(seed=None):
        if seed is None:
            scipy.stats.old_seed()
            return

        try:
            if len(seed) >= 2:
                scipy.stats.old_seed(seed[0], seed[1])
            else:
                scipy.stats.old_seed(seed[0])
        except TypeError:
            scipy.stats.old_seed(seed, 0)
    scipy.random = Bunch(seed = new_seed)

    # Move around fft function
    scipy.old_fft = scipy.fft
    scipy.fft = Bunch(rfft = scipy.old_fft, irfft = scipy.ifft)

    # Method renamed
    scipy.outer = scipy.outerproduct

    scipy.old_sum = scipy.sum
    def new_sum(x, axis=None, dtype=None, out=None):
        if (dtype is not None) or (out is not None):
            raise ValueError, "This use of sum is incompatible with old scipy."
        if (axis is not None):
            return scipy.old_sum(x, axis)

        output = x
        for ax in range(len(x.shape)):
            output = scipy.old_sum(output)
        return output
    scipy.sum = new_sum
