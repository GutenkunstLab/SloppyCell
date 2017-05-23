from SloppyCell.ReactionNetworks import *

from Nets import *

#
# The papers don't show any plots of the dynamics. The results below have been
#  checked against Ingenue as have all the underlying species not shown.
#

x0 = 0.1
thresh = lambda x: (x/x0)**3/(1+(x/x0)**3)

for ii, (net, t) in enumerate(zip(networks, int_times)):
    traj = Dynamics.integrate(net, [0, 1000])
    net_num = int(net.get_id()[4:])

    Plotting.figure(figsize=(10, 8))
    Plotting.clf()
    ticks = [0.05, 0.1, 0.2]
    Plotting.subplot(3, 1, 1)
    for spec_id in ['en_cell_%i' % jj for jj in range(4)]:
        Plotting.plot(traj.get_times(), thresh(traj.get_var_traj(spec_id)),
                      label = spec_id)
    Plotting.axis([0, 1000, 0, 1])
    Plotting.legend(loc='center right')
    Plotting.title('Stretched %s: score %f' % (net.get_id(), 
                                               param_scores[net_num]))
    Plotting.yticks(map(thresh, ticks), map(str, ticks))
    Plotting.hlines([thresh(0.1)], 0, 1000, 'k--')

    Plotting.subplot(3, 1, 2)
    for spec_id in ['wg_cell_%i' % jj for jj in range(4)]:
        Plotting.plot(traj.get_times(), thresh(traj.get_var_traj(spec_id)),
                      label = spec_id)
    Plotting.axis([0, 1000, 0, 1])
    Plotting.legend(loc='center right')
    Plotting.yticks(map(thresh, ticks), map(str, ticks))
    Plotting.hlines([thresh(0.1)], 0, 1000, 'k--')

    Plotting.subplot(3, 1, 3)
    for spec_id in ['hh_cell_%i' % jj for jj in range(4)]:
        Plotting.plot(traj.get_times(), thresh(traj.get_var_traj(spec_id)),
                      label = spec_id)
    Plotting.axis([0, 1000, 0, 1])
    Plotting.legend(loc='center right')
    Plotting.yticks(map(thresh, ticks), map(str, ticks))
    Plotting.hlines([thresh(0.1)], 0, 1000, 'k--')
