from SloppyCell.ReactionNetworks import *

from Nets import *

traj = Dynamics.integrate(net, [0, 100])
Plotting.figure(3)
Plotting.plot_trajectory(traj, ['C', 'M', 'X'])
