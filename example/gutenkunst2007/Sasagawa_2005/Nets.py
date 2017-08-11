from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000049.xml', 'base')
net.compile()

# Add some useful things for plotting
net.add_parameter('total_active_Rap1')
to_sum = [k for k in net.species.keys() if k.count('Rap1_GTP')]
net.add_assignment_rule('total_active_Rap1', ' + '.join(to_sum))

net.add_parameter('total_active_Ras')
to_sum = [k for k in net.species.keys() if k.count('Ras_GTP')]
net.add_assignment_rule('total_active_Ras', ' + '.join(to_sum))

net.add_parameter('total_active_ERK')
to_sum = [k for k in net.species.keys() if k.count('ppERK')]
net.add_assignment_rule('total_active_ERK', ' + '.join(to_sum))

# SBML file is for 10 ng/ml EGF
EGF_ng_ml = net.get_var_ic('EGF')/10.
net.add_parameter('EGF_ng_ml', EGF_ng_ml, is_optimizable=False)
# Convert using molecular weights
NGF_ng_ml = EGF_ng_ml * 6.1/26.
net.add_parameter('NGF_ng_ml', NGF_ng_ml, is_optimizable=False)

# Fixed 10 ng/ml EGF stimulation
EGF10_net = net.copy('EGF10')

# Ramp from 0 to 1.5 ng/ml EGF over 60 min
EGF1pt5_60_ramp = EGF10_net.copy('EGF1pt5_60_ramp')
EGF1pt5_60_ramp.add_assignment_rule('EGF', '1.5 * EGF_ng_ml * time/3600.')

# Fixed 10 ng/ml NGF stimulation
NGF10_net = net.copy('NGF10')
NGF10_net.set_var_ic('EGF', 0)
NGF10_net.set_var_ic('NGF', 10. * NGF_ng_ml)

# Ramp from 0 to 1.5 ng/ml NGF over 60 min
NGF1pt5_60_ramp = NGF10_net.copy('NGF1pt5_60_ramp')
NGF1pt5_60_ramp.add_assignment_rule('NGF', '1.5 * NGF_ng_ml * time/3600.')

networks = [EGF10_net, EGF1pt5_60_ramp, NGF10_net, NGF1pt5_60_ramp]
int_times = [(0, 60)] * 4
