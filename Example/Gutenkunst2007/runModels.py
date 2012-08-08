import os

import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *
from SloppyCell.ReactionNetworks.RunInParallel import *

def run_for_model(model):
    prev_dir = os.getcwd()
    os.chdir(model)
    sys.path.append(os.getcwd())
    import Nets

    # Update the typical values of all variables in each network to be the
    # maximum seen during integration of all the networks.
    PerfectData.update_typical_vals(Nets.networks, Nets.int_times)

    hessians = []
    key_sets = []

    for net, times in zip(Nets.networks, Nets.int_times):
        print 'Running for network %s.' % net.get_id()

        # A sensitivity with respect to log X doesn't make sens if X == 0, so
        # we don't calculate for those cases.
        for id in net.optimizableVars.keys():
            if net.get_var_ic(id) == 0:
                net.set_var_optimizable(id, False)

        # This model involve parameters that are exponents.
        # The resulting sensitivity equations have terms involving
        # log(species), which fails for species <= 0. This can cause problems
        # with the integration. To avoid those, we integrate with respect to
        # log(species). We also set all the initial conditions to be very
        # slighly greater than 0. One can check that the resulting sensitivity
        # trajectories still match those obtained from finite-difference on the
        # un-changed Network.
        # This model is a little finicky to integrate sensitivity for. It does fail on some
        # machines and not others.
        if model == 'von_Dassow_2000':
            net.integrateWithLogs = True
            for id in net.dynamicVars.keys():
                if net.get_var_ic(id) <= 0:
                    net.set_var_ic(id, scipy.misc.limits.double_resolution)

        if model in ['Chen_2004', 'von_Dassow_2000']:
            rtol = 1e-6
        else:
            rtol = 1e-9

        # We reset ics just in case
        net.resetDynamicVariables()
        sens_traj = Dynamics.integrate_sensitivity(net, times, rtol=rtol,
                                                   fill_traj=True)

        # This is specific to these two models, since they consider logarithmic
        # time.
        if model in ['Edelstein_1996', 'Brodersen_1987']:
            sens_traj.timepoints = scipy.log10(sens_traj.get_times())

        # We have 'data' on all species that aren't assigned constants.
        data_ids = [id for id in net.species.keys()
                    if not net.get_var_constant(id)]

        h, h_d = PerfectData.hessian_log_params(sens_traj, fixed_sf=True,
                                                return_dict=True,
                                                data_ids=data_ids)

        keys = sens_traj.optimizableVarKeys
        Utility.save((h, h_d, keys), 'hess_dict.%s.bp' % net.get_id())

        hessians.append(h)
        key_sets.append(keys)

        # Delete the trajectory we generated to save memory.
        del sens_traj
        del net.trajectory

    h, tot_keys = Utility.combine_hessians(hessians, key_sets)

    # Write hessian to file.
    scipy.savetxt('hessian.dat', h)
    f = file('hessian_keys.dat', 'w')
    f.write(os.linesep.join(tot_keys))
    f.close()

    sys.path.pop()
    del sys.modules['Nets']
    os.chdir(prev_dir)

import glob
import sys
if __name__ == '__main__':
    if len(sys.argv) != 1:
        directories = sys.argv[1:]
    else:
        import Common
        directories = [item[0] for item in Common.model_list]

    for dir in directories:
        print 'Running for model %s.' % dir
        run_for_model(dir)
