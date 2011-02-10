from SloppyCell.ReactionNetworks import *

from Nets import *

traj = Dynamics.integrate(prpp, [0, 70])
Plotting.subplot(2, 1, 1)
Plotting.plot_trajectory(traj, ['IMP'])
Plotting.axis([0, 70, 90, 120])
Plotting.subplot(2, 1, 2)
Plotting.plot_trajectory(traj, ['HX'])
Plotting.axis([0, 70, 5, 13])
