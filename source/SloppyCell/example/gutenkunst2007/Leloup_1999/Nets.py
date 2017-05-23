from SloppyCell.ReactionNetworks import *

net1 = IO.from_SBML_file('BIOMD0000000021.xml')
net1.add_assignment_rule('Pt', 'P0 + P1 + P2 + CC + Cn')
net1.add_assignment_rule('Tt', 'T0 + T1 + T2 + CC + Cn')

traj1 = Dynamics.integrate(net1, [0, 72, 84])
vals1 = traj1.get_var_vals(72)
net1.set_var_ics(vals1)

vals2 = traj1.get_var_vals(84)
net2 = net1.copy('net2')
net2.set_var_ics(vals2)

net3 = net1.copy('net3')
net3.set_var_ics(vals1)
net3.set_var_ic('V_mT', 0.28)
net3.set_var_ic('V_dT', 4.8)
traj3 = Dynamics.integrate(net3, [0, 200])
net3.set_var_ics(traj3.get_var_vals(200))

net4 = net1.copy('net4')
net4.set_var_ic('V_dT', 2.0)
net4.set_var_ic('V_mT', 0.99)
traj4 = Dynamics.integrate(net4, [0, 24*20 - 2])
net4.set_var_ics(traj4.get_var_vals_index(-1))

net5 = net4.copy('net5')
net4.set_var_ic('V_dT', 2.0)
net4.set_var_ic('V_mT', 0.99)
net5.set_var_ic('Mp', 1.8)
net5.set_var_ic('Mt', 0.5)
traj5 = Dynamics.integrate(net5, [0, 24*20 + 12])
net5.set_var_ics(traj5.get_var_vals_index(-1))

net6 = net1.copy('net6')
net6.set_var_ic('V_dT', 3.8)
net6.set_var_ic('V_mT', 0.4)
traj6 = Dynamics.integrate(net6, [0, 24*20 - 9])
net6.set_var_ics(traj6.get_var_vals_index(-1))

net7 = net6.copy('net7')
net7.set_var_ic('V_dT', 3.8)
net7.set_var_ic('V_mT', 0.4)
net7.set_var_ic('Mp', 0.2)
net7.set_var_ic('Mt', 4.0)
net7.set_var_ic('Cn', 1.0)
net7.set_var_ic('T2', 8)
traj7 = Dynamics.integrate(net7, [0, 24*20+24+6])
net7.set_var_ics(traj7.get_var_vals_index(-1))

networks = [net1, net2, net3, net4, net5, net6, net7]
int_times = [(0, 72), (0, 48), (0, 1000), (0, 96), (0, 96), (0, 96), (0, 96)]
