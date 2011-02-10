from SloppyCell.ReactionNetworks import *

from Nets import *

traj4a = Dynamics.integrate(net_12L_12D_12L_D, [0, 96])
traj6a = Dynamics.integrate(net_12D_L, [0, 60])

traj4b = Dynamics.integrate(mutant_12L_12D_12L_D, [0, 96])
traj6b = Dynamics.integrate(mutant_12D_L, [0, 60])

# Figure 4a
Plotting.figure(figsize=(8,10))
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(traj4a, ['cTm', 'cLm', 'light'])
Plotting.axis([0, 96, 0, 2.4])
Plotting.title('Figure 4, warning: scales in paper odd')
Plotting.subplot(2,1,2)
Plotting.plot_trajectory(traj4b, ['cTm', 'cLm', 'light'])
Plotting.axis([0, 96, 0, 7])

# Figure 6a
Plotting.figure(figsize=(8,10))
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(traj6a, ['cYm', 'light'])
Plotting.axis([0, 60, 0, 0.2])
Plotting.title('Figure 6, warning: scales in paper odd')
Plotting.subplot(2,1,2)
Plotting.plot_trajectory(traj6b, ['cYm', 'light'])
Plotting.axis([0, 60, 0, 0.6])
