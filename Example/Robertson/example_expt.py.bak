import numpy
from SloppyCell.ReactionNetworks import *
import example_net

traj_full = Dynamics.integrate(example_net.net, [0,40], fill_traj=True)

sample_traj = Dynamics.integrate(example_net.net, numpy.linspace(0,40,6), 
                                 fill_traj=False)

sigma = 0.25

numpy.random.seed(290137)

times = sample_traj.get_times()[1:]
data = {}
for var in ['A', 'B_scaled', 'C']:
    maxy = traj_full.get_var_traj(var).max()

    yy = sample_traj.get_var_traj(var)[1:]
    noisy_yy = yy*(1+sigma*numpy.random.randn(len(times)))
    uncert = sigma * maxy

    var_data = dict((t, (y, uncert)) for (t,y) in zip(times, noisy_yy))
    data[var] = var_data

expt = Experiment('example_expt',
                  {'example': data}, 
                  fixedScaleFactors={'A':1,'B_scaled':1, 'C':1})
