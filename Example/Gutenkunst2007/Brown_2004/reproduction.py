from SloppyCell.ReactionNetworks import *

from Nets import *

# Note that the plots in the paper have been rescaled, and that they refer to
#  an ensemble of parameters, not just the best fit parameters shown here.
# Thus the agreement is not exact.

# Integrate each of our networks
egf_stim_traj = Dynamics.integrate(egf_stim, [0, 120])
ngf_stim_traj = Dynamics.integrate(ngf_stim, [0, 120])
egf_ly_traj = Dynamics.integrate(egf_ly, [0, 120])
ngf_ly_traj = Dynamics.integrate(ngf_ly, [0, 120])
ngf_DN_Rap1_traj = Dynamics.integrate(ngf_DN_Rap1, [0, 120])
ngf_DN_Ras_traj = Dynamics.integrate(ngf_DN_Ras, [0, 120])

Plotting.figure(1)
Plotting.plot(egf_stim_traj.get_times(), 
              egf_stim_traj.get_var_traj('ErkActive'), 'r')
Plotting.plot(ngf_stim_traj.get_times(), 
              ngf_stim_traj.get_var_traj('ErkActive'), 'b')

Plotting.figure(2)
Plotting.subplot(2,1,1)
Plotting.plot(egf_ly_traj.get_times(),
              egf_ly_traj.get_var_traj('ErkActive'), 'r')
Plotting.plot(ngf_ly_traj.get_times(),
              ngf_ly_traj.get_var_traj('ErkActive'), 'b')

Plotting.subplot(2,1,2)
Plotting.plot(ngf_DN_Rap1_traj.get_times(),
              ngf_DN_Rap1_traj.get_var_traj('ErkActive'), 'r')
Plotting.plot(ngf_DN_Ras_traj.get_times(),
              ngf_DN_Ras_traj.get_var_traj('ErkActive'), 'b')
