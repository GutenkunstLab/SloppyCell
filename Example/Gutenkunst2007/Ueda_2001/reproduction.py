from __future__ import division
from past.utils import old_div
from SloppyCell.ReactionNetworks import *

from Nets import *

traj1 = Dynamics.integrate(net1, [0, 72])
Plotting.figure(1)
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(traj1, ['Perm', 'Perc', 'PTc', 'PTn'])
Plotting.axis([0, 72, 0, 6])
Plotting.subplot(2,1,2)
Plotting.plot_trajectory(traj1, ['Clkm', 'Clkc', 'CCn', 'CCc'])
Plotting.axis([0, 72, 0, 8])

traj2 = Dynamics.integrate(per01, [0, 72])
traj3 = Dynamics.integrate(dClkJrk, [0, 72])
traj4 = Dynamics.integrate(double, [0, 72])

maxP = max(traj1.get_var_traj('Perm'))
maxC = max(traj1.get_var_traj('Clkm'))
Plotting.figure(2)
Plotting.subplot(2,2,1)
Plotting.plot(traj1.get_times(), old_div(traj1.get_var_traj('Perm'),maxP))
Plotting.plot(traj1.get_times(), old_div(traj1.get_var_traj('Clkm'),maxC))
Plotting.axis([0, 72, 0, 1.5])
Plotting.subplot(2,2,2)
Plotting.plot(traj2.get_times(), old_div(traj2.get_var_traj('Perm'),maxP))
Plotting.plot(traj2.get_times(), old_div(traj2.get_var_traj('Clkm'),maxC))
Plotting.axis([0, 72, 0, 1.5])
Plotting.subplot(2,2,3)
Plotting.plot(traj3.get_times(), old_div(traj3.get_var_traj('Perm'),maxP))
Plotting.plot(traj3.get_times(), old_div(traj3.get_var_traj('Clkm'),maxC))
Plotting.axis([0, 72, 0, 1.5])
Plotting.subplot(2,2,4)
Plotting.plot(traj4.get_times(), old_div(traj4.get_var_traj('Perm'),maxP))
Plotting.plot(traj4.get_times(), old_div(traj4.get_var_traj('Clkm'),maxC))
Plotting.axis([0, 72, 0, 1.5])

traj5 = Dynamics.integrate(cper, [0, 72])
Plotting.figure(3)
Plotting.plot_trajectory(traj5, ['Perc', 'PTc', 'PTn'])
Plotting.axis([0,72, 0, 10])
