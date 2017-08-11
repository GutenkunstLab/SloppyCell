import os

import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *

import Common

Plotting.figure(figsize=(3.1,4))

for model_ii, (model, temp, temp) in enumerate(Common.model_list):
    print model
    # Load the hessian
    h = scipy.io.read_array(os.path.join(model, 'hessian.dat'))
    k = [item.strip() for item
         in file(os.path.join(model, 'hessian_keys.dat')).readlines()]

    e, v = Utility.eig(h)

    ii = -1
    while e[ii] == 0:
        non_zero_indices = scipy.compress(v[:,ii], range(len(v)))
        elems = [k[jj] for jj in non_zero_indices]
        print ('0 eigenvalue at index %i. '
               'Non-zero components of vector are %s' % (ii, str(elems)))
        ii -= 1

    # This takes care of any zero rows in h.  We count the corresponding 
    # parameter as P/I = 1.
    indices = []
    default_ones = 0
    for ii in range(h.shape[0]):
        if not scipy.all(h[ii] == 0):
            indices.append(ii)
        else:
            default_ones += 1

    pruned_h = h[indices,:]
    pruned_h = pruned_h[:,indices]
    P = scipy.diag(scipy.linalg.inv(pruned_h))
    I = 1/scipy.diag(pruned_h)
    yy = scipy.sort(scipy.sqrt(abs(I/P)))
    yy = scipy.concatenate((scipy.ones(default_ones), yy))
    # In once case, numerical noise gives us an apparent I/P > 1, which is
    # non-sensical. We manually make the fix here.
    yy = scipy.minimum(1, yy)

    width = (234/len(e))**0.25 * 0.25
    l = Plotting.plot_eigval_spectrum(yy, offset = 0.15+model_ii, 
                                      widths=0.7, lw=width)

ax = Plotting.gca()

ax.set_xticks(0.5 + scipy.arange(model_ii+1))
import string
xlabels = ['%s' % letter for letter in string.ascii_lowercase[:model_ii+1]]
ax.set_xticklabels(xlabels, fontsize=10, verticalalignment='bottom',
                   rotation=0, horizontalalignment='center')
for l in ax.get_xticklabels():
    l.set_y(l.get_position()[1] - 0.04)
for l in ax.get_xticklines():
    l.set_visible(False)
ax.set_xlim(0, model_ii+1)

ax.set_yscale('log',subsy=[])
ticks = range(-3, 1)
ax.set_yticks([10**ii for ii in ticks])
import matplotlib.lines
for l in ax.get_yticklines():
    if l.get_xdata() == (1,):
        l.set_visible(False)
    l.set_marker(matplotlib.lines.TICKLEFT)
ax.set_yticklabels([r'$10^{%s}$' % ii for ii in ticks], fontsize=10)
for l in ax.get_yticklabels():
    l.set_x(l.get_position()[0] - 0.04)
for ii in range(1, model_ii+1):
    ax.plot([ii, ii], [0.5e-3, 2e0], '-', lw=0.5, color=(0.75,0.75,0.75))
ax.set_xlim(0, model_ii+1)
ax.set_ylim(0.5e-3, 2e0)
Plotting.subplots_adjust(bottom=0.08, right=0.97, top=0.57)

Plotting.savefig('Fig1C_IoverP.eps')
