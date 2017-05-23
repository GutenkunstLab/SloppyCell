from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000051.xml', 'base')
net.compile()

traj = Dynamics.integrate(net, [0, 40])

Plotting.figure(1, figsize=(8,8))
Plotting.clf()
vars = ['cglcex','cfdp','cg1p','cg6p','cpep','cpyr','cf6p','cgap','cpg']
for ii, v in enumerate(vars):
    ax = Plotting.subplot(3,3, ii + 1)
    Plotting.plot_trajectory(traj, [v])
    if ax.is_last_row():
        Plotting.xticks([0, 10, 20, 30])
    else:
        Plotting.xticks([])
    if ii//3 == 0:
        Plotting.axis([0, 40, 0, 3.4])
    elif ii//3 == 1:
        Plotting.axis([0, 40, 0, 6.75])
    elif ii//3 == 2:
        Plotting.axis([0, 40, 0, 1.4])
    if not ax.is_first_col():
        Plotting.yticks([])
Plotting.subplots_adjust(wspace=0, hspace=0)

Plotting.figure(2, figsize=(6,10))
Plotting.clf()
vars = ['cg6p', 'cf6p', 'cpep', 'cpyr']
for ii, v in enumerate(vars):
    ax = Plotting.subplot(4,1, ii + 1)
    Plotting.plot_trajectory(traj, [v])
    if not ax.is_last_row():
        Plotting.xticks([])

    if ii == 0:
        Plotting.axis([0, 1, 3, 6])
    elif ii == 1:
        Plotting.axis([0, 1, 0.4, 0.9])
    elif ii == 2:
        Plotting.axis([0, 1, 1.5, 3.25])
    elif ii == 3:
        Plotting.axis([0, 1, 2.5, 4.25])
Plotting.subplots_adjust(wspace=0, hspace=0)

Plotting.figure(3, figsize=(6,10))
Plotting.clf()
vars = [('catp', 'cadp', 'camp'), ('cnad', 'cnadh'), ('cnadp', 'cnadph')]
for ii, v in enumerate(vars):
    ax = Plotting.subplot(3,1, ii + 1)
    Plotting.plot_trajectory(traj, v)
    if not ax.is_last_row():
        Plotting.xticks([])

    if ii == 0:
        Plotting.axis([0, 40, 0, 5])
    elif ii == 1:
        Plotting.axis([0, 40, 0, 2.0])
    elif ii == 2:
        Plotting.axis([0, 40, 0.04, 0.24])
Plotting.subplots_adjust(wspace=0, hspace=0)

