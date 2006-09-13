import scipy

if int(scipy.__version__.split('.')[1]) < 4:
    scipy.float_ = scipy.Float
    scipy.int_ = scipy.Int

    # Create an empty class to move around the limits module
    # http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/52308
    class Bunch(dict):
        def __init__(self,**kw):
            dict.__init__(self,kw)
            self.__dict__.update(kw)
    scipy.misc = Bunch(limits = scipy.limits)

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

    scipy.old_fft = scipy.fft
    scipy.fft = Bunch(rfft = scipy.old_fft, irfft = scipy.ifft)
