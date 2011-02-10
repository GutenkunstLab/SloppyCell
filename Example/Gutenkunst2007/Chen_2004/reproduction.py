from SloppyCell.ReactionNetworks import *

from Nets import *

traj = Dynamics.integrate(net, [0, 225])
traj_1 = Dynamics.integrate(mutant_1, [0, 225])
traj_2 = Dynamics.integrate(mutant_2, [0, 225])
traj_3 = Dynamics.integrate(mutant_3, [0, 225])
traj_4 = Dynamics.integrate(mutant_4, [0, 225])
traj_5 = Dynamics.integrate(mutant_5, [0, 225])
traj_6 = Dynamics.integrate(mutant_6, [0, 225])
traj_7 = Dynamics.integrate(mutant_7, [0, 350])
traj_8 = Dynamics.integrate(mutant_8, [0, 350])

Plotting.figure(1, figsize=(8,10.5))
Plotting.subplot(5, 1, 1)
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('BUD_3'), 'k:')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('ORI_43'), 'k--')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('SPN_53'), 'k-.')
Plotting.axis([0, max(traj.timepoints), 0, 3.2])

Plotting.subplot(5, 1, 2)
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CLB5T_23'), 'k:')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('PDS1_44'), 'k--')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CLB2T_21'), 'k-')
Plotting.axis([0, max(traj.timepoints), 0, 1.7])

Plotting.subplot(5, 1, 3)
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CDC6T_16'), 'k--')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('SIC1T_52'), 'k:')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CDH1_17'), 'k-')
Plotting.axis([0, max(traj.timepoints), 0, 1.4])

Plotting.subplot(5, 1, 4)
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CLN2_24'), 'k-')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CDC20_12'), 'k--')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('CDC14_8'), 'k:')
Plotting.axis([0, max(traj.timepoints), 0, 3])

Plotting.subplot(5, 1, 5)
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('SBF_49'), 'k:')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('MCM1_38'), 'k-')
Plotting.plot(traj.timepoints, traj.getVariableTrajectory('SWI5_54'), 'k--')
Plotting.axis([0, max(traj.timepoints), 0, 1.5])

Plotting.figure(2, figsize=(8,10.5))
Plotting.subplot(4, 2, 1)
Plotting.plot(traj_1.timepoints, traj_1.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_1.timepoints, traj_1.getVariableTrajectory('ORI_43'), 'k-.')
Plotting.plot(traj_1.timepoints, traj_1.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_1.timepoints, traj_1.getVariableTrajectory('CDH1_17'), 'k:')
Plotting.plot(traj_1.timepoints, traj_1.getVariableTrajectory('CDC6_14'), 'k-.')
#Plotting.title('cln1Del cln2Del cln3Del')
Plotting.axis([0, max(traj_1.timepoints), 0, 6])

Plotting.subplot(4, 2, 2)
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('ORI_43'), 'k-.')
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('CLB2T_21'), 'k--')
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('CDH1_17'), 'k:')
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('CDC6_14'), 'k-.')
Plotting.plot(traj_2.timepoints, traj_2.getVariableTrajectory('CLB5T_23'), 'k-')
#Plotting.title('cln1Del cln2Del cln3Del sic1Del')
Plotting.axis([0, max(traj_2.timepoints), 0, 6])

Plotting.subplot(4, 2, 3)
Plotting.plot(traj_3.timepoints, traj_3.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_3.timepoints, traj_3.getVariableTrajectory('ORI_43'), 'k-.')
Plotting.plot(traj_3.timepoints, traj_3.getVariableTrajectory('CLB2T_21'), 'k--')
Plotting.plot(traj_3.timepoints, traj_3.getVariableTrajectory('CDC6_14'), 'k-.')
Plotting.plot(traj_3.timepoints, traj_3.getVariableTrajectory('SIC1T_52'), 'k--')
#Plotting.title('cln1Del cln2Del cln3Del cdh1Del')
Plotting.axis([0, max(traj_3.timepoints), 0, 6])

Plotting.subplot(4, 2, 4)
Plotting.plot(traj_4.timepoints, traj_4.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_4.timepoints, traj_4.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_4.timepoints, traj_4.getVariableTrajectory('CDH1_17'), 'k:')
Plotting.plot(traj_4.timepoints, traj_4.getVariableTrajectory('CLB5T_23'), 
              'k--')
Plotting.plot(traj_4.timepoints, traj_4.getVariableTrajectory('CDC14_8'), 'k-.')
#Plotting.title('cdc20Del pds1Del')
Plotting.axis([0, max(traj_4.timepoints), 0, 4])

Plotting.subplot(4, 2, 5)
Plotting.plot(traj_5.timepoints, traj_5.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_5.timepoints, traj_5.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_5.timepoints, traj_5.getVariableTrajectory('CDH1_17'), 'k:')
Plotting.plot(traj_5.timepoints, traj_5.getVariableTrajectory('PDS1_44'), 'k-.')
Plotting.plot(traj_5.timepoints, traj_5.getVariableTrajectory('CDC14_8'), 'k-.')
#Plotting.title('cdc20Del clb5Del')
Plotting.axis([0, max(traj_5.timepoints), 0, 4])

Plotting.subplot(4, 2, 6)
Plotting.plot(traj_6.timepoints, traj_6.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_6.timepoints, traj_6.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_6.timepoints, traj_6.getVariableTrajectory('CDH1_17'), 'k:')
Plotting.plot(traj_6.timepoints, traj_6.getVariableTrajectory('CDC14_8'), 'k-.')
#Plotting.title('cdc20Del pds1Del clb5Del')
Plotting.axis([0, max(traj_6.timepoints), 0, 4])

Plotting.subplot(4, 2, 7)
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('CDC14_8'), 'k:')
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('CDC20_12'), 
              'k-.')
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('CDC14_8'), 'k-.')
Plotting.plot(traj_7.timepoints, traj_7.getVariableTrajectory('CLB5T_23'), 'k-')
#Plotting.title('sic1Del cdh1Del cdc6Del 2-49')
Plotting.axis([0, max(traj_8.timepoints), 0, 3])

Plotting.subplot(4, 2, 8)
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('MASS_37'), 'k-')
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('CLB2T_21'), 
              'k--')
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('CDC14_8'), 'k:')
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('CDC20_12'), 
              'k-.')
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('CDC14_8'), 'k-.')
Plotting.plot(traj_8.timepoints, traj_8.getVariableTrajectory('CLB5T_23'), 'k-')
#Plotting.title('sic1Del cdh1Del cdc6Del 2-49 GALL-CDC20')
Plotting.axis([0, max(traj_8.timepoints), 0, 3])
