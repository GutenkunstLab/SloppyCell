import scipy

from SloppyCell.ReactionNetworks import *

from Nets import *

specs = ['MF', 'MB', 'MD', 'MJ', 'MA', 'MC', 'ME', 'MH', 'MG', 'MK']
for var in specs:#net1.species.keys():
    net1.set_var_ic(var, net.get_var_ic(var))
traj = Dynamics.integrate(net1, [0, 80*60])

Plotting.figure()
for spec_id in specs:
    m0 = net.get_var_ic(spec_id)
    if spec_id == 'MH':
        # It looks very much like the variation in MH was multiplied by 50 for
        #  plotting purposes
        yy = 50*scipy.log(traj.get_var_traj(spec_id)/m0)
    else:
        # Also, MF and MD seem to be compared against 1, which isn't their
        #  steady-state deterministic value.
        yy = scipy.log(scipy.maximum(traj.get_var_traj(spec_id), 1)/max(m0, 1))
    Plotting.plot(traj.get_times()/60., yy, label=spec_id)

Plotting.legend()
Plotting.axis([0, 80, -2.5, 5])
