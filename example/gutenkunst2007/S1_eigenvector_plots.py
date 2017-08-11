import os
import string

import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *

import Common

for model_ii, (model, temp, temp) in enumerate(Common.model_list):
    # Load the hessian
    h = scipy.io.read_array(os.path.join(model, 'hessian.dat'))

    # Load the list of keys
    f = file(os.path.join(model, 'hessian_keys.dat'), 'r')
    keys = f.readlines()
    f.close()
    # Strip off extraneous characters
    keys = [k.strip() for k in keys]

    e, v = Utility.eig(h)
    v = scipy.real_if_close(v)

    Plotting.figure(figsize=(6, 4.5))
    N = 4
    # Plot the eigenvectors
    for ii in range(N):
        Plotting.subplot(N,1,ii+1)
        Plotting.plot_eigvect(v[:,ii], labels=keys)
        for p in Plotting.gca().patches:
            p.set_fc([0.5]*3)
        Plotting.gca().set_ylim(-1,1)
        Plotting.gca().set_yticks([-1,0,1])
        if ii == 0:
            Plotting.title('(%s) %s' % (string.ascii_lowercase[model_ii],model),
                           fontsize='x-large')

    # Strip the space out of the model name
    model_name = model.replace(' ', '')
    Plotting.savefig(os.path.join('eig_plots', '%s_eig.eps' % model_name))
    Plotting.close()
