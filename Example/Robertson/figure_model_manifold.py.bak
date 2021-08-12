"""
Plot Fig. 5
"""
from numpy import *
from SloppyCell.ReactionNetworks import *
import example_net, example_expt
import example_model

#
# Configure matplotlib fonts
#
import matplotlib
matplotlib.rc('font',**{'family':'sans-serif',
                        'sans-serif':['Arial'],
                        'style':'normal',
                        'size':10 })

sample_times = linspace(0,40,6)

#
# Calculate a model output surface
#
Npt = 41
land_xx = logspace(-4, 2, Npt)
land_yy = logspace(1.5, 7.5, Npt)
sample_data = empty((len(land_yy), len(land_xx), 3*(len(sample_times)-1)))
for ii,y in enumerate(land_yy):
    for jj,x in enumerate(land_xx):
        p = [x, y]

        # For late times, it helps to integrate with logs, which prevents
        # issues with the B trajectory crashing past 0. However,
        # we can't integrate with logs initiall, because several species
        # start at zero. So we break the interval up at delt.
        delt = 0.1
        traj = Dynamics.integrate(example_net.net, [0,delt], params=p, 
                                  fill_traj=False)
        # Continue integration                          
        newtimes = sample_times.copy()
        newtimes[0] = delt
        example_net.net.integrateWithLogs = True
        traj = Dynamics.integrate(example_net.net, newtimes, params=p,
                                  fill_traj=False)
        # Set back...
        example_net.net.integrateWithLogs = False

        # Stack trajectories of three interesting variables into single
        # row of data matrix.
        for kk, var in enumerate(['A','B_scaled','C']):
            start = (len(sample_times)-1)*kk
            end = (len(sample_times)-1)*(kk+1)
            sample_data[ii,jj,start:end] = abs(traj.get_var_traj(var)[1:])

# Take PCA in log-space of sample_data
# Note that eigenvectors are *columns* of v. And the direction with the 
# largest extent are the last entries.
reshaped = sample_data.reshape((-1, sample_data.shape[-1]))
# Exclude integrations that yielded a nan
reshaped = reshaped[logical_not(any(isnan(reshaped), axis=-1))]
u,v = Ensembles.PCA_eig(reshaped)
savetxt('example.model_manifold.nolog.PCA.us.dat', u)

# Project trajectories onto PCA axes
mean_sample_data = mean(reshaped, axis=0)
all_projs = dot(sample_data - mean_sample_data, v)

# Exclude nans
reshaped_all_projs = all_projs.reshape((-1, all_projs.shape[-1]))
reshaped_all_projs = reshaped_all_projs[logical_not(any(isnan(reshaped_all_projs), axis=-1))]

# Calculate projection of optimal trajectory
traj = Dynamics.integrate(example_net.net, sample_times, 
                          params=example_model.popt,
                          fill_traj=False)
popt_sample_data = zeros(all_projs.shape[-1])
for kk, var in enumerate(['A','B_scaled','C']):
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    popt_sample_data[start:end] = abs(traj.get_var_traj(var)[1:])
popt_proj = dot(popt_sample_data - mean_sample_data, v)

#
# Calculate projections of geodesics onto our manifold
#

geodesics = Utility.load('example.geodesics.bpkl')
all_gprojs = []
for xs,vs,ts in geodesics:
    gsamps = []
    for logp in xs:
        p = exp(logp)
        if not (1e-4 <= p[0] <= 1e2) or not (10**1.5 <= p[1] <= 10**7.5):
            continue

        # For late times, it helps to integrate with logs, which prevents
        # issues with the B trajectory crashing past 0. However,
        # we can't integrate with logs initiall, because several species
        # start at zero. So we break the interval up at delt.
        delt = 0.1
        traj = Dynamics.integrate(example_net.net, [0,delt], params=p, 
                                  fill_traj=False)
        # Continue integration                          
        newtimes = sample_times.copy()
        newtimes[0] = delt
        example_net.net.integrateWithLogs = True
        traj = Dynamics.integrate(example_net.net, newtimes, params=p,
                                  fill_traj=False)
        # Set back...
        example_net.net.integrateWithLogs = False

        # Stack trajectories of three interesting variables into single
        # row of data matrix.
        sample_data = []
        for kk, var in enumerate(['A','B_scaled','C']):
            start = (len(sample_times)-1)*kk
            end = (len(sample_times)-1)*(kk+1)
            sample_data.extend(abs(traj.get_var_traj(var)[1:]))
        gsamps.append(sample_data)
    gsamps = array(gsamps)
    gprojs = dot(gsamps - mean_sample_data, v)
    all_gprojs.append(gprojs)
        

#
# Calculate projection of data
#
data_sample_data = zeros(all_projs.shape[-1])
for kk, var in enumerate(['A','B_scaled','C']):
    data = [val for (time, (val, sigma)) in 
            sorted(example_expt.expt.get_data()['example'][var].items())]
    start = (len(sample_times)-1)*kk
    end = (len(sample_times)-1)*(kk+1)
    popt_sample_data[start:end] = data
data_proj = dot(popt_sample_data - mean_sample_data, v)

#
# Plot figure
#

import matplotlib.pyplot as pyplot
from mpl_toolkits.mplot3d import Axes3D

fig = pyplot.figure(290, figsize=(4.5,4.5), dpi=150)
fig.clear()
#
# Show PCAs, to illustrate how trajectories are changing.
#
for ii in range(1,4):
    ax = pyplot.subplot2grid((4,6),(0,2*(ii-1)), colspan=2, zorder=100)
    for jj, (var, color) in enumerate([('A','r'),('B','g'),('C','b')]):
        tt = sample_times[1:]
        l = len(tt)
        # Need to unpack sample row
        ax.plot(tt, v[jj*l:(jj+1)*l,-ii], '-o', ms=4, color=color)
    ax.set_xlabel('time')
    ax.set_xticks([0,20,40])
    ax.set_xlim(0,50)
    ax.set_title('PCA {0}'.format(ii), fontsize=10)
    ax.set_ylim(-1,1)
    ax.set_yticks([-1,-0.5,0,0.5,1])
    if ax.is_first_col():
        ax.set_ylabel(r'$\Delta y$')
    else:
        ax.set_yticklabels([])

ticks = [(-1,0,1,2),
         (-0.6, -0.3,0, 0.3),
         (-0.2,-0.1,0,0.1)]

## 3D plot of top three PCAs
ax3d = pyplot.subplot2grid((4,3), (1,0), colspan=2, rowspan=3, projection='3d',)
for ii in range(all_projs.shape[0]):
    ax3d.plot(all_projs[ii,:,-3], all_projs[ii,:,-2], 
              all_projs[ii,:,-1], '-k', linewidth=0.5, color=(0.25,0.25,0.25))
# Now plot lines connecting along x direction
for jj in range(all_projs.shape[1]):
    ax3d.plot(all_projs[:,jj,-3], all_projs[:,jj,-2], 
              all_projs[:,jj,-1], '-k', linewidth=0.5, color=(0.25,0.25,0.25))

#
# It's difficult to make the surface occlude things behind it. This attempts
# to do it with patches, but it's doesn't really work well.
#
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#verts = []
#for ii in range(all_projs.shape[0])[:-1]:
#    for jj in range(all_projs.shape[1])[:-1]:
#        verts.append([all_projs[ii,jj,-3:], 
#                      all_projs[ii,jj+1,-3:], 
#                      all_projs[ii+1,jj,-3:], 
#                      ])
#        verts.append([all_projs[ii+1,jj+1,-3:], 
#                      all_projs[ii,jj+1,-3:], 
#                      all_projs[ii+1,jj,-3:], 
#                      ])
#poly = Poly3DCollection(verts, facecolors='w', edgecolors='w', linewidth=0.5)
#ax3d.add_collection(poly)

# Plot edges more boldly
ax3d.plot(all_projs[:,-1,-3], all_projs[:,-1,-2], 
          all_projs[:,-1,-1], '-k', linewidth=2)
ax3d.plot(all_projs[:,0,-3], all_projs[:,0,-2], 
          all_projs[:,0,-1], '-k', linewidth=2)
ax3d.plot(all_projs[0,:,-3], all_projs[0,:,-2], 
          all_projs[0,:,-1], '-k', linewidth=2)
ax3d.plot(all_projs[-1,:,-3], all_projs[-1,:,-2], 
          all_projs[-1,:,-1], '-k', linewidth=2)
ax3d.plot([popt_proj[-3]], [popt_proj[-2]], [popt_proj[-1]], 'or', ms=6, mec='k', zorder=1000)
ax3d.plot([data_proj[-3]], [data_proj[-2]], [data_proj[-1]], '*b', ms=8, zorder=2000)
ax3d.view_init(39, 64)
ax3d.set_zlabel('PCA 1')
ax3d.set_zticks(ticks[0])
ax3d.set_zlim(ticks[0][0], ticks[0][-1])
ax3d.set_ylabel('PCA 2')
ax3d.set_yticks(ticks[1])
ax3d.set_ylim(ticks[1][0], ticks[1][-1])
ax3d.set_xlabel('PCA 3')
ax3d.set_xticks(ticks[2])
ax3d.set_xlim(xmin=ticks[2][0], xmax=0.14)

for gprojs in all_gprojs:
    ax3d.plot(gprojs[:,-3], gprojs[:,-2], gprojs[:,-1], '-r', lw=1)

# 2D plots of top PCAs, for clarity
for ii in range(1,4):
    ax = pyplot.subplot2grid((4,6), (ii,4), colspan=2)
    # Trick to get pairs of 1-2, 1-3, 2-3
    indices = range(1,4)
    del indices[3-ii]
    for ii in range(all_projs.shape[0]):
        ax.plot(all_projs[ii,:,-indices[0]], all_projs[ii,:,-indices[1]], 
                '-k', linewidth=0.5, color=(0.25,0.25,0.25))
    # Now plot lines connecting along x direction
    for jj in range(all_projs.shape[1]):
        ax.plot(all_projs[:,jj,-indices[0]], all_projs[:,jj,-indices[1]], 
                '-k', linewidth=0.5, color=(0.25,0.25,0.25))
    ax.plot(all_projs[:,-1,-indices[0]], all_projs[:,-1,-indices[1]], '-k', 
              linewidth=1)
    ax.plot(all_projs[:,0,-indices[0]], all_projs[:,0,-indices[1]], '-k', 
              linewidth=1)
    ax.plot(all_projs[-1,:,-indices[0]], all_projs[-1,:,-indices[1]], '-k', 
              linewidth=1)
    ax.plot(all_projs[0,:,-indices[0]], all_projs[0,:,-indices[1]], '-k', 
              linewidth=1)
    ax.plot([popt_proj[-indices[0]]], [popt_proj[-indices[1]]], 'or', ms=4, zorder=100, mec='k')
    ax.plot([data_proj[-indices[0]]], [data_proj[-indices[1]]], '*b', ms=6, zorder=2000)
    for gprojs in all_gprojs:
        ax.plot(gprojs[:,-indices[0]], gprojs[:,-indices[1]], '-r', lw=0.75)
    ax.set_xlabel('PCA {0}'.format(indices[0]))
    ax.set_ylabel('PCA {0}'.format(indices[1]))
    ax.set_xticks(ticks[indices[0]-1])
    ax.set_yticks(ticks[indices[1]-1])

fig.tight_layout()
# tight_layout isn't doing great here
# Let's tweak additionally
fig.subplots_adjust(wspace=0.8, hspace=0.73, top=0.95)
pos = ax3d.get_position()
ax3d.set_position([0.01,0.05,pos.x0+pos.width*0.97,pos.y0+pos.height*1.2])
ax3d.set_zorder(-10)

fig.text(0.055,0.94,'A')
fig.text(0.405,0.94,'B')
fig.text(0.69,0.94,'C')
fig.text(0.025,0.71, 'D')
fig.text(0.645,0.71, 'E')
fig.text(0.645,0.47, 'F')
fig.text(0.645,0.23, 'G')

tr = fig.transFigure
r = matplotlib.patches.Rectangle((0,0.74),1,0.26,fill=True,color='w',transform=fig.transFigure)
fig.patches.append(r)

fig.savefig('nonlog_model_manifold.pdf')

pyplot.show()
