from SloppyCell.ReactionNetworks import *

from Nets import *

traj1 = Dynamics.integrate(net1, [0, 150*60])
Plotting.figure(1)
Plotting.plot_trajectory(traj1, ['MAPK_PP', 'MAPK'])

traj2 = Dynamics.integrate(net2, [0, 205*60])
Plotting.figure(2)
Plotting.plot_trajectory(traj2, ['MAPK_PP', 'MAPK'])
