from SloppyCell.ReactionNetworks import *

from Nets import *

A_traj = Dynamics.integrate(A, [0, 20])

Plotting.figure(1)
Plotting.plot_trajectory(A_traj, ['L'], loc='lower right')

B_traj = Dynamics.integrate(B, [0, 20])
C_traj = Dynamics.integrate(C, [0, 20])
D_traj = Dynamics.integrate(D, [0, 20])

Plotting.figure(1)
Plotting.subplot(2,2,1)
Plotting.plot_trajectory(A_traj, ['L'], loc='lower right')
Plotting.axis([0, 20, 0, 1])
Plotting.subplot(2,2,2)
Plotting.plot_trajectory(B_traj, ['L'])
Plotting.axis([0, 20, 0, 1])
Plotting.subplot(2,2,3)
Plotting.plot_trajectory(C_traj, ['L_p'], loc='lower right')
Plotting.axis([0, 20, 0, 1])
Plotting.subplot(2,2,4)
Plotting.plot_trajectory(D_traj, ['L_p'], loc='upper right')
Plotting.axis([0, 20, 0, 1])

E_traj = Dynamics.integrate(E, [0, 20])
F_traj = Dynamics.integrate(F, [0, 40])
G_traj = Dynamics.integrate(G, [0, 20])
H_traj = Dynamics.integrate(H, [0, 20])

Plotting.figure(2)
Plotting.subplot(2,2,1)
Plotting.plot_trajectory(E_traj, ['D'])
Plotting.axis([0, 20, 0, 1])
Plotting.subplot(2,2,2)
Plotting.plot_trajectory(F_traj, ['D'])
Plotting.axis([0, 40, 0, 1])
Plotting.subplot(2,2,3)
Plotting.plot_trajectory(G_traj, ['W'])
Plotting.axis([0, 20, 0, 1])
Plotting.subplot(2,2,4)
Plotting.plot_trajectory(H_traj, ['W'])
Plotting.axis([0, 20, 0, 1])
