#!/usr/bin/env python
#PBS -l nodes=1,pmem=2048mb,mem=4096mb,ncpus=1,cput=6:00:00
#PBS -m abe
#PBS -M rng7@cornell.edu

# For batch running
import os
if os.getenv('PBS_O_WORKDIR'):
    os.chdir(os.getenv('PBS_O_WORKDIR'))
    os.sys.path.append(os.getcwd())

from SloppyCell.ReactionNetworks import *

from LeeNet import net

net.set_var_ic('W', 0)
traj1 = Dynamics.integrate(net, [0, 1e5], rtol=1e-12)
net1 = net.copy('unstimulated')
net1.set_var_ics(traj1.get_var_vals_index(-1))

net.set_var_ic('W', 1)
traj2 = Dynamics.integrate(net, [0, 1e5], rtol=1e-12)
net2 = net.copy('stimulated')
net2.set_var_ics(traj2.get_var_vals_index(-1))

net3 = net1.copy('transient_a')
net3.add_parameter('lam', 1./20, is_optimizable=False)
net3.add_assignment_rule('W', 'exp(-lam*time)')
traj3 = Dynamics.integrate(net3, [0, 16*60])

net4 = net3.copy('transient_b')
net4.set_var_ic('v14', net3.get_var_ic('v14')*5)
net4.set_var_ic('k15', net3.get_var_ic('k15')*5)
traj4 = Dynamics.integrate(net4, [0, 16*60])

net5 = net3.copy('transient_c')
net5.set_var_ic('v14', net3.get_var_ic('v14')/5)
net5.set_var_ic('k15', net3.get_var_ic('k15')/5)
traj5 = Dynamics.integrate(net5, [0, 16*60])

fig2a = net1.copy('fig2a')
fig2a.set_var_ic('v12', 0)
fig2a.set_var_optimizable('v12', False)
fig2a.set_var_ic('v14', 0)
fig2a.set_var_optimizable('v14', False)
traj2a = Dynamics.integrate(fig2a, [0, 3*60], rtol=1e-12)

fig2b = net1.copy('fig2b')
fig2b.set_var_ic('v12', 0)
fig2b.set_var_optimizable('v12', False)
fig2b.set_var_ic('v14', 0)
fig2b.set_var_optimizable('v14', False)
fig2b.set_var_ic('X12', 0.2)
traj2b = Dynamics.integrate(fig2b, [0, 3*60], rtol=1e-12)

fig2c = net1.copy('fig2c')
fig2c.set_var_ic('v12', 0)
fig2c.set_var_optimizable('v12', False)
fig2c.set_var_ic('v14', 0)
fig2c.set_var_optimizable('v14', False)
fig2c.set_var_ic('k2', 0)
fig2c.set_var_optimizable('k2', False)
fig2c.set_var_ic('X2', 1000)
fig2c.set_var_ic('Dsh0', 1100)
traj2c = Dynamics.integrate(fig2c, [0, 3*60], rtol=1e-12)

fig2d = net1.copy('fig2d')
fig2d.set_var_ic('v12', 0)
fig2d.set_var_optimizable('v12', False)
fig2d.set_var_ic('v14', 0)
fig2d.set_var_optimizable('v14', False)
fig2d.set_var_ic('k4', 0)
fig2d.set_var_optimizable('k4', False)
fig2d.set_var_ic('k9', 0)
fig2d.set_var_optimizable('k9', False)
traj2d = Dynamics.integrate(fig2d, [0, 3*60], rtol=1e-12)

###
### This doesn't appear to be correct...
###
fig2e = net1.copy('fig2e')
fig2e.set_var_ic('v12', 0)
fig2e.set_var_optimizable('v12', False)
fig2e.set_var_ic('v14', 0)
fig2e.set_var_optimizable('v14', False)
fig2e.set_var_ic('TCF0', 2015)
traj2e = Dynamics.integrate(fig2e, [0, 3*60], rtol=1e-12)

for traj in [traj2a, traj2b, traj2c, traj2d, traj2e]:
    Plotting.plot(traj.get_times(), traj.get_var_traj('BCatenin'))


import pd
networks = [net3, net4, net5]
int_times = [[0, 16*60]] * 3
#pd.update_typical_vals(networks, int_times)
#pd.run_for_nets2(networks, int_times, send_email=True, rtol=1e-9)

#Plotting.figure(1, (5, 10))
#Plotting.subplot(2,1,1)
#Plotting.plot(traj3.get_times()/60., traj3.get_var_traj('BCatenin'))
#Plotting.plot(traj4.get_times()/60., traj4.get_var_traj('BCatenin'))
#Plotting.plot(traj5.get_times()/60., traj5.get_var_traj('BCatenin'))
#Plotting.axis([0, 16, 34, 72])
#Plotting.subplot(2,1,2)
#Plotting.plot(traj3.get_times()/60., 1000*traj3.get_var_traj('Axin'))
#Plotting.plot(traj4.get_times()/60., 1000*traj4.get_var_traj('Axin'))
#Plotting.plot(traj5.get_times()/60., 1000*traj5.get_var_traj('Axin'))
#Plotting.axis([0, 16, 0, 22])

#import compare_sens_mod as cm
##sens_traj = Utility.load('sens.transient_a.bp')
#for id in net.optimizableVars.keys():
#    if id != 'GSK0':
#        net3.set_var_optimizable(id, False)
#
#sens_traj = Dynamics.integrate_sensitivity_2(net3, [0, 1000])
#cm.compare_sens(5, 0, net3, sens_traj)
