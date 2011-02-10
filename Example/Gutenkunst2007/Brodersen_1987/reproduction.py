import scipy

from SloppyCell.ReactionNetworks import *

from Nets import *

alb_net.resetDynamicVariables()
alb_times = scipy.logspace(-6, -2, 1000)
alb_traj = Dynamics.integrate(alb_net, alb_times) 

heme_net.resetDynamicVariables()
heme_times = scipy.logspace(-1, 3, 1000)
heme_traj = Dynamics.integrate(heme_net, heme_times)

Plotting.figure(1, figsize=(8,10))
Plotting.subplot(2,1,1)
Plotting.plot_trajectory(alb_traj, ['ligand_per_protein'], logx=True)
Plotting.axis([1e-6, 1e-2, 0, 6])
Plotting.title('Albumin best fit')

Plotting.subplot(2,1,2)
Plotting.plot_trajectory(heme_traj, ['ligand_per_protein'], logx=True)
Plotting.axis([1e-1, 1e3, 0, 4])
Plotting.title('Hemoglobin best fit')
