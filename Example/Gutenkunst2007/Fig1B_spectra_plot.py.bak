import os

import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *

import Common

Plotting.figure(figsize=(3.1,4))

for model_ii, (model, temp, temp) in enumerate(Common.model_list):
    # Load the hessian
    h = scipy.io.read_array(os.path.join(model, 'hessian.dat'))

    e, v = Utility.eig(h)
    e = scipy.real(e)

    width = (234/len(e))**0.25 * 0.25
    l = Plotting.plot_eigval_spectrum(e/max(e), offset = 0.15+model_ii, 
                                      widths=0.7, lw=width)

# Now a lot of fiddling to make the plot prettier
ax = Plotting.gca()
for ii in range(1, model_ii+1):
    ax.plot([ii, ii], [0.5e-6, 2], '-', lw=0.5, color=(0.75,0.75,0.75))

# Add labels
ax.set_xticks(0.5 + scipy.arange(model_ii+1))
import string
xlabels = ['%s' % tup for tup in string.ascii_lowercase[:model_ii+1]]
                                             
ax.set_xticklabels(xlabels, fontsize=10, verticalalignment='bottom',
                   rotation=0, horizontalalignment='center')
for l in ax.get_xticklabels():
    l.set_y(l.get_position()[1] - 0.04)
for l in ax.get_xticklines():
    l.set_visible(False)
ax.set_xlim(0, model_ii+1)

ax.set_yscale('log',subsy=[])
ticks = range(-6, 1)
ax.set_yticks([10**ii for ii in ticks])
import matplotlib.lines
for l in ax.get_yticklines():
    if l.get_xdata() == (1,):
        l.set_visible(False)
    l.set_marker(matplotlib.lines.TICKLEFT)
ax.set_yticklabels([r'$10^{%s}$' % ii for ii in ticks], fontsize=10)
for l in ax.get_yticklabels():
    l.set_x(l.get_position()[0] - 0.04)
ax.set_ylim(0.5e-6, 2)
Plotting.subplots_adjust(bottom=0.08, right=0.97, top=0.97)

Plotting.savefig('Fig1B_spectra_plot.eps')
