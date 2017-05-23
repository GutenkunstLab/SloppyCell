import scipy
import scipy.io

from SloppyCell.ReactionNetworks import *

h = scipy.io.read_array('Brown_2004/hessian.dat')
h_rs = scipy.io.read_array('Brown_2004_Rescaled/hessian.dat')
h_a1 = scipy.io.read_array('Brown_2004_All_One/hessian.dat')

Plotting.figure(figsize=(4,3))
u, v = Utility.eig(h)
#l = Plotting.plot_eigvals(u/max(u), join=True)
Plotting.plot_eigval_spectrum(u/max(u), lw=1)
u, v = Utility.eig(h_rs)
#l_rs = Plotting.plot_eigvals(u/max(u), join=True)
Plotting.plot_eigval_spectrum(u/max(u), offset=1.2, lw=1)
u, v = Utility.eig(h_a1)
Plotting.plot_eigval_spectrum(u/max(u), offset=2.4, lw=1)
#l_a1 = Plotting.plot_eigvals(u/max(u), join=True)

Plotting.gca().set_xlim(-.2, 3.6)
Plotting.gca().set_ylim(0.5e-6, 2)
#Plotting.legend((l, l_rs, l_a1), ('Original model', 'Rescaled', 'All One'))
Plotting.subplots_adjust(left=0.16, bottom=0.16)
Plotting.xticks([0.5, 1.7, 2.9], ['Original', 'Rescaled', 'All One'], fontsize='medium')
Plotting.ylabel(r'$\lambda/\lambda_0$', fontsize='large')
Plotting.savefig('S4_rescaled_Brown_plot.eps')
