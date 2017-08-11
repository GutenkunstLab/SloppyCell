from scipy import *
import scipy.io as io

from SloppyCell.ReactionNetworks import *

alb_dat = io.read_array('Brodersen_1987/Albumin.ens.dat')
alb_dat = transpose(reshape(alb_dat, (6, 30)))
ua,va = Ensembles.PCA_eig_log_params(alb_dat)

alb_chi2_hess, temp, temp = Utility.load('Brodersen_1987/hess_dict.albumin.bp')
ua2, va2 = Utility.eig(alb_chi2_hess)

heme_dat = io.read_array('Brodersen_1987/Hemoglobin.ens.dat')
heme_dat = transpose(reshape(heme_dat, (4, 30)))
uh,vh = Ensembles.PCA_eig_log_params(heme_dat)

heme_chi2_hess, temp, temp = Utility.load('Brodersen_1987/hess_dict.hemeglobin.bp')
uh2, vh2 = Utility.eig(heme_chi2_hess)

Plotting.figure(figsize=(6,3))
Plotting.plot_eigval_spectrum(ua/ua[0], widths=0.8)
Plotting.plot_eigval_spectrum(ua2/ua2[0], widths=0.8, offset=1.0)
Plotting.plot_eigval_spectrum(uh/uh[0], widths=0.8, offset=2.0)
Plotting.plot_eigval_spectrum(uh2/uh2[0], widths=0.8, offset=3.0)
Plotting.gca().set_xlim(-0.1, 3.9)
Plotting.gca().set_ylim(0.5e-6, 2)
Plotting.ylabel(r'$\lambda/\lambda_0$', fontsize='large')
Plotting.xticks([0.4, 1.4, 2.4, 3.4],
                ('Alb ens', 'Alb chi^2', 'Heme ens', 'Heme chi^2'))

Plotting.savefig('S5_Brodersen.eps')
