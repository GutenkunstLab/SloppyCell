from SloppyCell.ReactionNetworks import *

from Nets import *

traj_EGF10 = Dynamics.integrate(EGF10_net, [0, 60*60])
traj_EGF1pt5_60_ramp = Dynamics.integrate(EGF1pt5_60_ramp, [0, 60*60])
traj_NGF10 = Dynamics.integrate(NGF10_net, [0, 60*60])
traj_NGF1pt5_60_ramp = Dynamics.integrate(NGF1pt5_60_ramp, [0, 60*60])

Plotting.figure(3)
Plotting.clf()
Plotting.plot(traj_EGF10.get_times()/60., 
              traj_EGF10.get_var_traj('total_active_Ras'))
Plotting.plot(traj_NGF10.get_times()/60., 
              traj_NGF10.get_var_traj('total_active_Ras'))
Plotting.axis([0, 60, 0, 0.1])

Plotting.figure(4)
Plotting.clf()
Plotting.plot(traj_EGF10.get_times()/60., 
              traj_EGF10.get_var_traj('total_active_Rap1'))
Plotting.plot(traj_NGF10.get_times()/60., 
              traj_NGF10.get_var_traj('total_active_Rap1'))
Plotting.axis([0, 60, 0, 0.042])

Plotting.figure(8)
Plotting.clf()
traj_EGF1pt5_60_ramp
Plotting.plot(traj_EGF1pt5_60_ramp.get_times()/60., 
              100*traj_EGF1pt5_60_ramp.get_var_traj('total_active_ERK')/net.get_var_ic('ERK'))
Plotting.axis([0, 60, 0, 100])

Plotting.figure(9)
Plotting.clf()
Plotting.plot(traj_NGF1pt5_60_ramp.get_times()/60., 
              100*traj_NGF1pt5_60_ramp.get_var_traj('total_active_ERK')/net.get_var_ic('ERK'))
Plotting.axis([0, 60, 0, 100])
