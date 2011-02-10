from SloppyCell.ReactionNetworks import *

import Nets
reload(Nets)
from Nets import *

base_traj = Dynamics.integrate(base, (0, 400))
r_traj = Dynamics.integrate(reduced, [0, 200])
stable_traj = Dynamics.integrate(stable, [0, 400])

Plotting.figure(2)
Plotting.subplot(2,1,1)
Plotting.plot(base_traj.get_times(), base_traj.get_var_traj('A'), '-k')
Plotting.axis([0, 400, 0, 2000])
Plotting.subplot(2,1,2)
Plotting.plot(base_traj.get_times(), base_traj.get_var_traj('R'), '-k')
Plotting.axis([0, 400, 0, 2000])

Plotting.figure(3)
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(r_traj, ['R', 'C'])
Plotting.axis([0, 200, 0, 2000])
Plotting.subplot(2,1,2)
Plotting.plot_trajectory(base_traj, ['R', 'C'])
Plotting.axis([0, 200, 0, 3000])

Plotting.figure(5)
Plotting.plot_trajectory(stable_traj, ['R'])
Plotting.axis([0, 400, 0, 3000])
