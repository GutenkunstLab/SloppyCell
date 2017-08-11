from SloppyCell.ReactionNetworks import *

net1 = IO.from_SBML_file('BIOMD0000000022.xml')
net1.set_var_constant('EmptySet', True)
traj = Dynamics.integrate(net1, [0, 20*24])
vals = traj.get_var_vals(20*24)
net1.set_var_ics(vals)

per01 = net1.copy('per01')
per01.set_var_ic('s2', 0)
traj = Dynamics.integrate(per01, [0, 200*24])
vals = traj.get_var_vals(200*24)
per01.set_var_ics(vals)

dClkJrk = net1.copy('dClkJrk')
dClkJrk.set_var_ic('s6', 0)
traj = Dynamics.integrate(dClkJrk, [0, 20*24])
vals = traj.get_var_vals(20*24)
dClkJrk.set_var_ics(vals)

double = net1.copy('double')
double.set_var_ic('s6', 0)
traj = Dynamics.integrate(double, [0, 20*24])
vals = traj.get_var_vals(20*24)
double.set_var_ics(vals)

cper = net1.copy('cper')
cper.set_var_ic('s1', 0)
cper.set_var_ic('c1', 0.846)
traj = Dynamics.integrate(cper, [0, 20*24+15])
vals = traj.get_var_vals(20*24+15)
cper.set_var_ics(vals)

# per01, dClkJrk, and double are all in fixed-points. Thus we neglect them,
#  since their dynamics are very uninteresting.
networks = [cper, net1]
int_times = [(0, 72)] * len(networks)
