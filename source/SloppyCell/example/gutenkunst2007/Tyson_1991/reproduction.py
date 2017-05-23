from SloppyCell.ReactionNetworks import *

from Nets import *

#####
## Figure 3a from the paper:
#####

# Make sure our networks look correct
traj1 = Dynamics.integrate(base_net, [0, 100])
Plotting.figure(1)
Plotting.plot_trajectory(traj1, ['M', 'YT'])
Plotting.axis([0, 100, 0, .4])

#####
## Figure 3b from the paper:
#####

traj2 = Dynamics.integrate(perturbed_net, [0, 100])
Plotting.figure(2)
Plotting.plot_trajectory(traj2, ['M', 'YT'], logy=True)
Plotting.axis([0, 100, 1e-4, 1e0])

#####
## Figure 3c from paper
#####

traj3 = Dynamics.integrate(growth_net, [0, 500])
Plotting.figure(4)
Plotting.plot_trajectory(traj3, ['M', 'YT'])
Plotting.axis([0, 500, 0, 0.6])

Plotting.figure(5)
Plotting.plot_trajectory(traj3, ['k6'])
Plotting.axis([0, 500, 1, 3])
