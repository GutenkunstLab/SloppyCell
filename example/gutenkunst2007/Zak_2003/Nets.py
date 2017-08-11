from SloppyCell.ReactionNetworks import *

# First we need to convert our model to level 2, version 1
IO.SBML.to_SBML_l2v1('Model.xml', 'Model-l2v1.xml')
net = IO.from_SBML_file('Model-l2v1.xml', 'Zak2003')

# I don't know why these are even in the Network.
net.remove_component('time')
net.remove_component('time_source')
net.remove_component('Rtime')

# Run to a steady state.
net.set_var_ic('S_t', 0)
traj = Dynamics.integrate(net, [0, 2000*60])
vals = traj.get_var_vals_index(-1)
net.set_var_ics(vals)

net1 = net.copy('single_pulse')
net1.add_event('begin_pulse', 'gt(time, 60*10.)', {'S_t': 1})
net1.add_event('end_pulse', 'gt(time, 60*10+10.)', {'S_t': 0})

net2 = net.copy('double_pulse')
net2.add_event('begin_pulse', 'gt(time, 60*10.)', {'S_t': 1})
net2.add_event('end_pulse', 'gt(time, 60*10+60.)', {'S_t': 0})
net2.add_event('begin_pulse_2', 'gt(time, 60*10. + 24*60)', {'S_t': 1})
net2.add_event('end_pulse_2', 'gt(time, 60*10+60. + 24*60)', {'S_t': 0})

net3 = net.copy('step')
net3.add_event('begin_step', 'gt(time, 60*10.)', {'S_t': 1})

networks = [net1, net2, net3]
int_times = [[0, 48*60]] * len(networks)
