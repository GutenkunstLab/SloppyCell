import os

import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *

import Common

pts_per_param = 100
# Fraction uncertainty assumed for data points
f = 0.1

Plotting.figure(figsize=(8/2.54, 13/2.54))
Plotting.gcf().set_facecolor('w')
Plotting.axes([0, 0, 1, 1])

for model_ii, (model, N_c, N_s) in enumerate(Common.model_list):
    print model
    h = scipy.io.read_array(os.path.join(model, 'hessian.dat'))

    e, v = Utility.eig(h)
    
    N_d = pts_per_param * h.shape[0]
    h *= N_d/(f**2 * N_c * N_s)

    # There is one case where a parameter has no effect in the quadratic
    #  approximation. To account for this, we count rows of the hessian that
    #  are all zeros, then filter them out.
    indices = []
    default_ones = 0
    for ii in range(h.shape[0]):
        if not scipy.all(h[ii] == 0):
            indices.append(ii)
        else:
            default_ones += 1

    pruned_h = h[indices,:]
    pruned_h = pruned_h[:,indices]
    h = pruned_h

    sigma_log = scipy.sqrt(abs(scipy.diag(scipy.linalg.pinv(h, 1e-19))))
    sigma_log = scipy.sort(sigma_log)
    interval_size = scipy.exp(4*sigma_log) - 1
    # We assign any parameter whose uncertainty we can't calculate an 
    #  uncertainty of zero.
    interval_size = scipy.concatenate((interval_size, 
                                       scipy.ones(default_ones) * scipy.inf))
    yy = scipy.where(interval_size > 1e5, 1e5, interval_size)

    bins = [0, .1, 1, 10, 100, scipy.inf]
    left_add = 0.1 + model_ii/9*(len(bins)+0.25)
    left = scipy.arange(len(bins) - 1) + left_add
    height = []
    for bottom, top in zip(bins[:-1], bins[1:]):
        height.append(scipy.sum((bottom < yy) & (yy < top))/(1.0*h.shape[0]))
    bottom = 10-(model_ii%9) * 1.4
    Plotting.bar(left, height, width=1.0, bottom = bottom)
    label = model + ' (%i)' % h.shape[0]
    Plotting.text(left_add + (len(bins)-1)*0.5, bottom-0.1, label, 
                  horizontalalignment='center',
                  verticalalignment='top', fontsize=10)

bottom = 10 - 7*1.4 - 0.5
labels = ('< 10%','10% - 1', '1 - 10', '10 - 100', '> 100')
for xx, label in zip(left, labels):
    Plotting.text(xx + 0.3, bottom - 0.05, label, verticalalignment='top', horizontalalignment='left', rotation=-30, fontsize=10)
Plotting.gca().get_frame().set_visible(False)
Plotting.gca().set_xticks([])
Plotting.gca().set_yticks([])

Plotting.savefig('Fig3.eps')
