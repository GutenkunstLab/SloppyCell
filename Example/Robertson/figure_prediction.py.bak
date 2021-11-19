"""
Plot Fig. 4
"""
import string
from SloppyCell.ReactionNetworks import *
from numpy import *
import example_model

import matplotlib.pyplot as pyplot

#
# Configure matplotlib fonts
#
import matplotlib
matplotlib.rc('font',**{'family':'sans-serif',
                        'sans-serif':['Arial'],
                        'style':'normal',
                        'size':10 })


def prune_ens(ens, Npt):
    # Prune data ensemble to Npt elements
    return ens[1::len(ens)//Npt]

# Max dimensions are 4.5 x 7.25 inches
pyplot.close(210)
fig = pyplot.figure(210, figsize=(4.5,2.0), dpi=150)
fig.clear()

(ens_data, ens_all_measured, ens_one_measured, ens_jtj) =\
    Utility.load('example.ensembles.bpkl')
(pred_data_lw, pred_data_up, 
 pred_all_measured_lw, pred_all_measured_up,
 pred_one_measured_lw, pred_one_measured_up, 
 pred_jtj_lw, pred_jtj_up) = Utility.load('example.pred_bounds.bpkl')
tt = linspace(0,50,100)
land_xx, land_yy, Z = Utility.load('example.model_surface.bpkl')
Z[isnan(Z)] = Z[logical_not(isnan(Z))].max()

ens_list = [ens_all_measured, ens_one_measured, ens_jtj, ens_data]
pred_list = [(pred_all_measured_lw, pred_all_measured_up),
             (pred_one_measured_lw, pred_one_measured_up),
             (pred_jtj_lw, pred_jtj_up),
             (pred_data_lw, pred_data_up)]
for ii, (ens, (pred_lw, pred_up)) in enumerate(zip(ens_list, pred_list)):
    ax_land = pyplot.subplot2grid((2,4),(0,ii), aspect=1)
    ax_pred = pyplot.subplot2grid((2,4), (1,ii))

    norm = matplotlib.colors.LogNorm(vmin=example_model.cost_opt, vmax=10**2.5)
    mappable = ax_land.pcolor(land_xx[::2], land_yy[::2], Z[::2,::2], norm=norm, cmap='gray')
    ax_land.set_xscale('log', subsx=[])
    ax_land.set_yscale('log', subsy=[])
    # Note that these limits are square
    ax_land.set_xlim(1e-4, 1e2)
    ax_land.set_ylim(10**1.5, 10**7.5)
    
    ax_land.set_yticks([1e2, 1e7])
    ax_land.set_ylabel('$k_3$', horizontalalignment='left')
    ax_land.set_xticks([1e-4,1e2])
    ax_land.set_xlabel('$k_1$', verticalalignment='bottom')

    ax_land.text(0.09, 0.91, string.uppercase[2*ii], ha='left', va='top', 
                 transform=ax_land.transAxes)
        

    # Add pruned data ensemble (100 samples)
    ens_pruned = prune_ens(ens, 100)
    ax_land.plot(ens_pruned[:,0], ens_pruned[:,1], 'ow', ms=3, mew=0.5)
    
    ax_pred.plot(tt, pred_up, '-k', lw=0.5)
    ax_pred.plot(tt, pred_lw, '-k', lw=0.5)
    ax_pred.fill_between(tt, pred_lw, pred_up, color=(0.5,0.5,0.5))
    ax_pred.set_xlabel('time')
    ax_pred.set_xticks([0, 20, 40])
    ax_pred.set_ylim(0,0.50)
    ax_pred.set_yticks([0,0.25,0.5])
    if ax_pred.is_first_col():
        ax_pred.set_ylabel('[C]')
    else:
        ax_pred.set_yticklabels([])
    ax_pred.text(0.03, 0.95, string.uppercase[2*ii+1], ha='left', va='top', 
                 transform=ax_pred.transAxes)
    ax_pred.set_xlim(0,50)

fig.tight_layout()
fig.subplots_adjust(bottom=0.19, left=.13, right=0.97, top=0.95, wspace=0.26, hspace=0.67)

pyplot.show()

fig.savefig('../figures/ensemble_predictions.pdf')
