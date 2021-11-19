from SloppyCell.ReactionNetworks import *

# Modifications to SBML...
# Removed function LD, because it used 'ceil' which is not something we can deal
#  with
# Replaced variable value_of_LD with light (more descriptive name)
# Replaced calls to LD with light
# Removed timeOfDay and dayLength variables
net = IO.from_SBML_file('BIOMD055-noceil.xml', 'base')
net.compile()

# Set up a network that will switch light on/off at 12 hour intervals.
net1212 = net.copy('net_1212')
net1212.set_var_ic('light', 1)
net1212.add_parameter('turntime', 12, is_constant=False)
net1212.add_event('light_switch', 'gt(time, turntime)', {'light': '1-light',
                                                         'turntime': '12+time'})
mutant_net = net1212.copy('cca1lhy')
mutant_net.set_var_ic('p1', net.get_var_ic('p1')/1000)

# Run to the limit cycle 
traj = Dynamics.integrate(net1212, [0, 24*10])
net1212.set_var_ics(traj.get_var_vals_index(-1))
# Go to limit cycle
traj = Dynamics.integrate(mutant_net, [0, 24*10])
mutant_net.set_var_ics(traj.get_var_vals_index(-1))

net_12L_12D_12L_D = net1212.copy('net_12L_12D_12L_D')
net_12L_12D_12L_D.remove_component('light_switch')
net_12L_12D_12L_D.remove_component('turntime')
net_12L_12D_12L_D.set_var_ic('light', 1)
net_12L_12D_12L_D.add_event('off_12', 'gt(time, 12)', {'light': 0})
net_12L_12D_12L_D.add_event('on_24', 'gt(time, 24)', {'light': 1})
net_12L_12D_12L_D.add_event('off_36', 'gt(time, 36)', {'light': 0})

# Run for twelve more hours to get to the dark part of the cycle
traj = Dynamics.integrate(net1212, [0, 12])
net1212.set_var_ics(traj.get_var_vals_index(-1))

net_12D_L = net1212.copy('net_12D_L')
net_12D_L.remove_component('light_switch')
net_12D_L.remove_component('turntime')
net_12D_L.set_var_ic('light', 0)
net_12D_L.add_event('on_12', 'gt(time, 12)', {'light': 1})

mutant_12L_12D_12L_D = mutant_net.copy('mutant_12L_12D_12L_D')
mutant_12L_12D_12L_D.remove_component('light_switch')
mutant_12L_12D_12L_D.remove_component('turntime')
mutant_12L_12D_12L_D.set_var_ic('light', 1)
mutant_12L_12D_12L_D.add_event('off_12', 'gt(time, 12)', {'light': 0})
mutant_12L_12D_12L_D.add_event('on_24', 'gt(time, 24)', {'light': 1})
mutant_12L_12D_12L_D.add_event('off_36', 'gt(time, 36)', {'light': 0})
trajm = Dynamics.integrate(mutant_12L_12D_12L_D, [0, 96])

# Run for twelve more hours to get to the dark part of the cycle
traj = Dynamics.integrate(mutant_net, [0, 12])
mutant_net.set_var_ics(traj.get_var_vals_index(-1))

mutant_12D_L = mutant_net.copy('mutant_12D_L')
mutant_12D_L.remove_component('light_switch')
mutant_12D_L.remove_component('turntime')
mutant_12D_L.set_var_ic('light', 0)
mutant_12D_L.add_event('on_12', 'gt(time, 12)', {'light': 1})

networks = [net_12L_12D_12L_D, net_12D_L, mutant_12L_12D_12L_D, mutant_12D_L]
int_times = [(0, 96), (0, 96), (0, 48), (0,48)]
