import scipy

from SloppyCell.ReactionNetworks import *

from Nets import *

net1.resetDynamicVariables()
traj1 = Dynamics.integrate(net1, (10**-5, 10**2))

net2.resetDynamicVariables()
traj2 = Dynamics.integrate(net2, (10**-3, 10**4))

Plotting.figure(1, figsize=(8,10))
Plotting.subplot(3,1,1)
Plotting.semilogx(traj1.get_times(), 
              1 - (traj1.get_var_traj('A') + traj1.get_var_traj('AL')
                   + traj1.get_var_traj('ALL')))
Plotting.axis([10**-5, 10**2, 0, 1])
Plotting.subplot(3,1,2)
Plotting.plot_trajectory(traj1, ['B', 'BL', 'ALL', 'ILL', 'DLL'], logx=True)
Plotting.axis([10**-5, 10**2, 0, 1])

Plotting.subplot(3,1,3)
Plotting.plot_trajectory(traj2, ['B', 'IL', 'ALL', 'ILL', 'DLL', 'DL', 
                                 'I', 'D'], logx=True)
Plotting.axis([10**-3, 10**4, 0, 1])

#
# We don't use the conditions in the following figures for the hessian 
#  calculations, since they're very repetative, and don't add much more
#  information than the progression and recovery networks.
#

# Figure 8: a and b
Plotting.figure(4, figsize=(8,10))
Plotting.subplot(2,1,1)
for Lconc in scipy.logspace(-6, -3, 16):
    net1.set_var_ic('L', Lconc)
    net1.resetDynamicVariables()
    traj1 = Dynamics.integrate(net1, (10**-6, 10**1))
    Plotting.semilogx(traj1.get_times(), 1 - (traj1.get_var_traj('A') + traj1.get_var_traj('AL') + traj1.get_var_traj('ALL')))

Plotting.subplot(2,1,2)
net3 = net1.copy('lowerK2')
net3.set_var_ic('kf_11', 200)
for Lconc in scipy.logspace(-6, -3, 16):
    net3.set_var_ic('L', Lconc)
    net3.resetDynamicVariables()
    traj1 = Dynamics.integrate(net3, (10**-6, 10**1))
    Plotting.semilogx(traj1.get_times(), 1 - (traj1.get_var_traj('A') + traj1.get_var_traj('AL') + traj1.get_var_traj('ALL')))

# Figure 9: a and b
Plotting.figure(6, figsize=(8, 10))
dnet = net.copy('4em7')
dnet.set_var_ic('L', 4e-7)
dnet.resetDynamicVariables()
traj = Dynamics.integrate(dnet, [10**-6, 10**3])
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(traj, logx=True)
Plotting.axis([10**-6, 10**3, 0, 1])
Plotting.subplot(2,1,2)
Plotting.plot_trajectory(traj)
Plotting.axis([-100, 1000, 0, 1])
