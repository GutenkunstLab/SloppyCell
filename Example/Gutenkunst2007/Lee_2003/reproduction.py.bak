from SloppyCell.ReactionNetworks import *

import Nets

traj2a = Dynamics.integrate(Nets.fig2a, [0, 3*60])
traj2b = Dynamics.integrate(Nets.fig2b, [0, 3*60])
traj2c = Dynamics.integrate(Nets.fig2c, [0, 3*60])
traj2d = Dynamics.integrate(Nets.fig2d, [0, 3*60])
traj2e = Dynamics.integrate(Nets.fig2e, [0, 3*60])

Plotting.figure(2)
for traj in [traj2a, traj2b, traj2c, traj2d, traj2e]:
    percent = 100*traj.get_var_traj('BCatenin')/traj.get_var_val('BCatenin', 0)
    Plotting.plot(traj.get_times()/60., percent, '-k')
    Plotting.axis([0, 3, 0, 105])

traj6a = Dynamics.integrate(Nets.fig6a, [0, 16*60])
traj6b = Dynamics.integrate(Nets.fig6b, [0, 16*60])
traj6c = Dynamics.integrate(Nets.fig6c, [0, 16*60])

Plotting.figure(6, (5, 10))
Plotting.subplot(2,1,1)
Plotting.plot(traj6a.get_times()/60., traj6a.get_var_traj('BCatenin'), '-k')
Plotting.plot(traj6b.get_times()/60., traj6b.get_var_traj('BCatenin'), '-r')
Plotting.plot(traj6c.get_times()/60., traj6c.get_var_traj('BCatenin'), '-g')
Plotting.axis([-1, 16, 34, 72])
Plotting.subplot(2,1,2)
Plotting.plot(traj6a.get_times()/60., 1000*traj6a.get_var_traj('Axin'), '-k')
Plotting.plot(traj6b.get_times()/60., 1000*traj6b.get_var_traj('Axin'), '-r')
Plotting.plot(traj6c.get_times()/60., 1000*traj6c.get_var_traj('Axin'), '-g')
Plotting.axis([-1, 16, 0, 22])
