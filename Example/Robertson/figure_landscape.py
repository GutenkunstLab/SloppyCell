"""
Plots Fig. 2
"""
from SloppyCell.ReactionNetworks import *
import scipy.stats
from numpy import *
import example_model

import matplotlib.pyplot as pyplot
pyplot.close(110)
#
# Configure matplotlib fonts
#
import matplotlib
matplotlib.rc('font',**{'family':'sans-serif',
                        'sans-serif':['Arial'],
                        'style':'normal',
                        'size':10 })

# Max dimensions are 4.5 x 7.25 inches
fig = pyplot.figure(110, figsize=(3.6,5.8), dpi=200)
fig.clear()

# 
# This is a complicated figure to create, because we want axes on a non-regular
# grid, and we want to ensure one is square. Thus we can't use the normal
# subplot-based layout.
#

# Give us the whole figure to place plots into.
fig.subplots_adjust(left=0,bottom=0,top=1,right=1)

# We use a gridspec to have complete control over plot placement. The grid is
# scaled so each unit is 1/100 of an inch.
gridheight = int(fig.get_figheight()*125)
gridwidth = int(fig.get_figwidth()*125)
gs = matplotlib.gridspec.GridSpec(gridheight, gridwidth)

# Create the axes we'll use. I did it all at once to make it easier to quickly
# lay out.
ax_traj = fig.add_subplot(gs[10:160, 100:-100])
# It's helpful to add labels initiall, because they'll affect the space neeeded
# between axes
ax_traj.set_xlabel('time')
# Top marginal distribution. 
ax_top = fig.add_subplot(gs[220:290,65:355])
ax_top.set_xticklabels([])
ax_top.set_yticks([])
# Central landscape plot. Note that I was careful to ensure it's square.
ax_land = fig.add_subplot(gs[300:590,65:355])
ax_land.set_xlabel('$k_1$', labelpad=0)
ax_land.set_ylabel('$k_3$', labelpad=0)
# Right marginal distribution
ax_rt = fig.add_subplot(gs[300:590,365:435])
ax_rt.set_yticklabels([])
ax_rt.set_xticks([])
# Colorbar axes
cax = fig.add_subplot(gs[650:665,125:305])

# Labels for figure panels
fig.text(0.12,0.96, 'A', fontsize=14)
fig.text(0.07,0.67, 'B', fontsize=14)

#
# Plot of model trajectory and data
#
# Note that colors are chosen to match Eydgahi plots
net = example_model.m.calcColl[0]
traj = Dynamics.integrate(net, [0, 45], params=example_model.popt)
for var, color, factor in [('A','r',1),('B_scaled','g',1),('C','b',1)]:
    ax_traj.plot(traj.get_times(), traj.get_var_traj(var)*factor, color)
    expt = example_model.m.exptColl.values()[0]
    tt = expt.data['example'][var].keys()
    yy = [expt.data['example'][var][t][0]*factor for t in tt]
    yerr = [expt.data['example'][var][t][1]*factor for t in tt]
    ax_traj.errorbar(tt, yy, yerr, linestyle='o', color='k')
    ax_traj.plot(tt, yy, 'o', color=color, ms=5)

# Label lines
ax_traj.text(1.3, 1.0, '[A]', color='r')
ax_traj.text(3.6, 0.45, '[B]', color='g')
ax_traj.text(2.0, 0.30, '(x 10 )', color='g')
# Horribly hacky way to put in a superscript.
ax_traj.text(9.0, 0.41, '4', color='g', fontsize=6, 
             horizontalalignment='center', verticalalignment='center')
ax_traj.text(32, 0.35, '[C]', color='b')
ax_traj.set_ylim(0, 1.3)
ax_traj.set_yticks([0, 0.5, 1.0])

#
# Landscape
#
land_xx, land_yy, Z = Utility.load('example.model_surface.bpkl')
# Fill in values where calculation failed.
Z[isnan(Z)] = Z[logical_not(isnan(Z))].max()
norm = matplotlib.colors.LogNorm(vmin=example_model.cost_opt, vmax=10**2.5)
ens_data = Utility.load('example.ens_data.bpkl')
mappable = ax_land.pcolor(land_xx, land_yy, Z, norm=norm, cmap='gray')
ax_land.set_xscale('log')
ax_land.set_yscale('log')
# Note that these limits are square
ax_land.set_xlim(1e-4, 1e2)
ax_land.set_ylim(10**1.5, 10**7.5)

def prune_ens(ens, Npt):
    # Prune data ensemble to Npt elements
    return ens[1::len(ens)//Npt]
# Add pruned data ensemble (100 samples)
ens_data_pruned = prune_ens(ens_data, 100)
ax_land.plot(ens_data_pruned[:,0], ens_data_pruned[:,1], 'ow', ms=4, zorder=10)

# Best-fit
ax_land.plot([example_model.popt[0]], [example_model.popt[1]], 'or', zorder=30)

def contour_top_bottom(h, delC=1, Npts=100):
    """
    Contours tracing ellipse from the hessian approximation

    Note that these are countours in delta X and delta Y.
    """
    delx = sqrt(delC*h[1,1]/(h[1,1]*h[0,0] - h[0,1]**2))
    xx = linspace(-delx, delx, Npts)

    a = h[1,1]
    b = 2*h[0,1]*xx
    c = h[0,0]*xx**2 - delC

    # Avoid errors due to tiny negative values
    inner = b**2 - 4*a*c
    inner = maximum(inner, 0)

    top = (-b + sqrt(inner))/(2*a)
    bottom = (-b - sqrt(inner))/(2*a)
    
    return xx, top, bottom

# Hessian 95% confidence interval
delxx,delyy_top,delyy_bottom = contour_top_bottom(example_model.jtj, 
                                                  delC=4, Npts=100)
# Convert from deltas in log params to normal params
xx = exp(log(example_model.popt[0]) + delxx)
yy_top = exp(log(example_model.popt[1]) + delyy_top)
yy_bottom = exp(log(example_model.popt[1]) + delyy_bottom)
# This concatenation ensures we plot a single smooth curve.
xx = concatenate((xx, xx[::-1], [xx[0]]))
yy = concatenate((yy_top, yy_bottom[::-1], [yy_top[0]]))
ax_land.plot(xx, yy, '-b', zorder=20)

geodesics = Utility.load('example.geodesics.bpkl')
for xs, vs, ts in geodesics:
    ax_land.plot(exp(xs[:,0]), exp(xs[:,1]), '-r', lw=1.5, zorder=2)

#
# Colorbar scale for landscape
#
cbar = fig.colorbar(mappable, cax=cax, orientation='horizontal')
cbar.set_ticks(concatenate([[example_model.cost_opt],10**arange(1,2.51,0.5)]))
cbar.set_ticklabels([r'${0:.1f}$'.format(example_model.cost_opt)] + 
                    [r'$10^{{{0:.1f}}}$'.format(_) for _ in arange(1,2.51,0.5)])
for tl in cbar.ax.get_xticklabels():
    tl.set_va('bottom')
    tl.set_y(-1.2)
cbar.set_label('cost', labelpad=0)

#
# Top marginal distribution
#
# Match extend of axes fo landscape plot
xmin, xmax = ax_land.get_xlim()
# Note that we bin in terms of log-params.
bins = linspace(log10(xmin), log10(xmax), 31, endpoint=False)
ax_top.hist(log10(ens_data[:,0]), bins=bins, normed=True, fc='w')
# pdf in terms of log-params
scale = example_model.jtj_uncerts[0]*log10(exp(1))
pdf = scipy.stats.distributions.norm(loc=log10(example_model.popt[0]), scale=scale).pdf
xx = linspace(log10(xmin), log10(xmax), 100)
yy = pdf(xx)
ax_top.plot(xx, yy, '-b')
ax_top.set_xlim(log10(xmin), log10(xmax))
ax_top.set_ylim(0, yy.max()*1.7)

#
# Right marginal distribution
#
# This is a little tricky because x and y are swapped.
ymin, ymax = ax_land.get_ylim()
bins = linspace(log10(ymin), log10(ymax), 31, endpoint=False)
q = histogram(log10(ens_data[:,1]), bins, normed=True)[0]
ax_rt.barh(bottom=bins[:-1], width=q, height=diff(bins), fc='w')
scale = example_model.jtj_uncerts[1]*log10(exp(1))
pdf = scipy.stats.distributions.norm(loc=log10(example_model.popt[1]), scale=scale).pdf
yy = linspace(log10(ymin), log10(ymax), 100)
xx = pdf(yy)
ax_rt.plot(xx, yy, '-b')
ax_rt.set_ylim(log10(ymin), log10(ymax))
ax_rt.set_xlim(0, xx.max()*1.7)

pyplot.show()

fig.savefig('../figures/ensemble_fits.pdf')
