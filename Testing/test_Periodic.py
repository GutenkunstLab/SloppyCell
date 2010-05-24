import unittest
import numpy

from SloppyCell.ReactionNetworks import *

# Example Drosophila circadian network from:
# Leloup & Goldbeter. Journal of Biological Rhythms 13: 70 (1999).

plotTraj = False # Recreate Figure 4 from the paper if True
if plotTraj:
    try:
        import pylab
        pylab.figure(figsize=(16, 11.425))
        pylab.subplots_adjust(top=0.95, bottom=0.05, hspace=0.35)
    except: plotTraj = False

net = Network('test')
net.add_compartment('cell')

parameters = {'v_sP': 0.8, 'v_sT': 1.0, 
              'v_mP': 0.8, 'v_mT': 0.7, 
              'K_mP': 0.2, 'K_mT': 0.2, 
              'k_sP': 0.9, 'k_sT': 0.9, 
              'v_dP': 2.0, 'v_dT': 2.0, 
              'k_1': 1.2, 'k_2': 0.2, 'k_3': 1.2, 'k_4': 0.6, 
              'K_IP': 1.0, 'K_IT': 1.0, 
              'K_dP': 0.2, 'K_dT': 0.2, 
              'n': 4, 
              'K_1P': 2.0, 'K_2P': 2.0, 'K_3P': 2.0, 'K_4P': 2.0,
              'K_1T': 2.0, 'K_2T': 2.0, 'K_3T': 2.0, 'K_4T': 2.0, 
              'k_d': 0.01, 'k_dC': 0.01, 'k_dN': 0.01, 
              'V_1P': 8.0, 'V_1T': 8.0, 
              'V_2P': 1.0, 'V_2T': 1.0, 
              'V_3P': 8.0, 'V_3T': 8.0, 
              'V_4P': 1.0, 'V_4T': 1.0}

for name, val in parameters.items():
    net.add_parameter(name, val)

net.add_species('M_P', 'cell')
net.add_species('P_0', 'cell')
net.add_species('P_1', 'cell')
net.add_species('P_2', 'cell')
net.add_species('M_T', 'cell')
net.add_species('T_0', 'cell')
net.add_species('T_1', 'cell')
net.add_species('T_2', 'cell')
net.add_species('C', 'cell')
net.add_species('C_N', 'cell')

# Transcription
net.addReaction('M_P transcription',
                kineticLaw='v_sP*(K_IP**n)/((K_IP**n)+(C_N)**n)',
                stoichiometry={'M_P':+1, 'C_N':0})
net.addReaction('M_T transcription',
                kineticLaw='v_sT*(K_IT**n)/((K_IT**n)+(C_N)**n)',
                stoichiometry={'M_T':+1, 'C_N':0})

# Translation
net.addReaction('P_0 translation', kineticLaw='k_sP*M_P',
                stoichiometry={'M_P':0, 'P_0':+1})
net.addReaction('T_0 translation', kineticLaw='k_sT*M_T',
                stoichiometry={'M_T':0, 'T_0':+1})

# Phosphorylation
net.addReaction('P_0 phosphorylation', kineticLaw='V_1P*P_0/(K_1P+P_0)',
                stoichiometry={'P_0':-1, 'P_1':+1})
net.addReaction('P_1 dephosphorylation', kineticLaw='V_2P*P_1/(K_2P+P_1)',
                stoichiometry={'P_1':-1, 'P_0':+1})
net.addReaction('P_1 phosphorylation', kineticLaw='V_3P*P_1/(K_3P+P_1)',
                stoichiometry={'P_1':-1, 'P_2':+1})
net.addReaction('P_2 dephosphorylation', kineticLaw='V_4P*P_2/(K_4P+P_2)',
                stoichiometry={'P_1':+1, 'P_2':-1})

net.addReaction('T_0 phosphorylation', kineticLaw='V_1T*T_0/(K_1T+T_0)',
                stoichiometry={'T_0':-1, 'T_1':+1})
net.addReaction('T_1 dephosphorylation', kineticLaw='V_2T*T_1/(K_2T+T_1)',
                stoichiometry={'T_1':-1, 'T_0':+1})
net.addReaction('T_1 phosphorylation', kineticLaw='V_3T*T_1/(K_3T+T_1)',
                stoichiometry={'T_1':-1, 'T_2':+1})
net.addReaction('T_2 dephosphorylation', kineticLaw='V_4T*T_2/(K_4T+T_2)',
                stoichiometry={'T_1':+1, 'T_2':-1})

# Complex formation and transport
net.addReaction('P-T binding', kineticLaw='k_3*P_2*T_2',
                stoichiometry={'P_2':-1, 'T_2':-1, 'C':+1})
net.addReaction('C unbinding', kineticLaw='k_4*C',
                stoichiometry={'P_2':+1, 'T_2':+1, 'C':-1})
net.addReaction('C nuclear transport', kineticLaw='k_1*C',
                stoichiometry={'C':-1, 'C_N':+1})
net.addReaction('C_N cytoplasmic transport', kineticLaw='k_2*C_N',
                stoichiometry={'C':+1, 'C_N':-1})

# Specific degradation
net.addReaction('M_P degradation', kineticLaw='v_mP*M_P/(K_mP+M_P)',
                stoichiometry={'M_P':-1})
net.addReaction('M_T degradation', kineticLaw='v_mT*M_T/(K_mT+M_T)',
                stoichiometry={'M_T':-1})

net.addReaction('P_2 degradation', kineticLaw='v_dP*P_2/(K_dP+P_2)',
                stoichiometry={'P_2':-1})
net.addReaction('T_2 degradation', kineticLaw='v_dT*T_2/(K_dT+T_2)',
                stoichiometry={'T_2':-1})

# Non-specific degradation
net.addReaction('M_P nonspecific degradation', kineticLaw='k_d*M_P',
                stoichiometry={'M_P':-1})
net.addReaction('P_0 nonspecific degradation', kineticLaw='k_d*P_0',
                stoichiometry={'P_0':-1})
net.addReaction('P_1 nonspecific degradation', kineticLaw='k_d*P_1',
                stoichiometry={'P_1':-1})
net.addReaction('P_2 nonspecific degradation', kineticLaw='k_d*P_2',
                stoichiometry={'P_2':-1})
net.addReaction('M_T nonspecific degradation', kineticLaw='k_d*M_T',
                stoichiometry={'M_T':-1})
net.addReaction('T_0 nonspecific degradation', kineticLaw='k_d*T_0',
                stoichiometry={'T_0':-1})
net.addReaction('T_1 nonspecific degradation', kineticLaw='k_d*T_1',
                stoichiometry={'T_1':-1})
net.addReaction('T_2 nonspecific degradation', kineticLaw='k_d*T_2',
                stoichiometry={'T_2':-1})

net.addReaction('C nonspecific degradation', kineticLaw='k_dC*C',
                stoichiometry={'C':-1})
net.addReaction('C_N nonspecific degradation', kineticLaw='k_dN*C_N',
                stoichiometry={'C_N':-1})

net.add_species('P_t', 'cell')
net.add_assignment_rule('P_t', 'P_0+P_1+P_2+C+C_N')
net.add_species('T_t', 'cell')
net.add_assignment_rule('T_t', 'T_0+T_1+T_2+C+C_N')

base_net = net.copy('test_Periodic')

base_net.compile()

DDic = {'M_P': 0.339861640394,  'P_0': 0.084951960872,
        'P_1': 0.0761208617382, 'P_2': 0.0379109639304,
        'M_T': 1.00941245929,   'T_0': 0.282665734452,
        'T_1': 0.261767573957,  'T_2': 0.161161625329,
        'C': 0.0959575363908,   'C_N': 0.799716886506}

LDic = {'M_P': 0.091103911665979753, 'P_0': 0.023248542312634367,
        'P_1': 0.022096471574764676, 'P_2': 0.01251106143540503,
        'M_T': 1.4252746160225265,   'T_0': 0.54130471189106055,
        'T_1': 0.799209240863616,    'T_2': 4.7317657065068524,
        'C': 0.17887971663488947,    'C_N': 1.202640586603108}

class test_Periodic(unittest.TestCase):
    def test_Unforced(self):
        """Test that an unforced limit cycle can be found with the correct phase"""
        period = 25.036999481467504

        net = base_net.copy('test_Unforces')

        for id, val in DDic.items(): net.set_var_ic(id, round(val, 2))

        net.set_periodic(period=round(period,2), phase=('T_t', 'min'),
                         xtol=0.00001, log=False, minVel=0.0001)

        net.Calculate({'M_P':[0, 72]})
        traj = net.trajectory

        if plotTraj:
            sp = pylab.subplot(321)
            for var in ['M_P', 'P_t', 'M_T', 'T_t']:
                pylab.plot(traj.timepoints, traj.get_var_traj(var),
                           label=var)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(0, 10)
            pylab.yticks([0,2,4,6,8,10])
            pylab.ylabel('Total PER and TIM proteins (P_t and T_t)\n'
                         'and their mRNAs (M_P and M_T)', fontsize='small',
                         horizontalalignment='center')
            pylab.legend(loc=0)
            pylab.title('A')

            sp = pylab.subplot(323)
            pylab.plot(traj.timepoints, traj.get_var_traj('M_P'), label='M_P')
            pylab.plot(traj.timepoints, (traj.get_var_traj('M_T')-0.5)/3.5*1.4,
                       label='M_T')
            pylab.plot(traj.timepoints, traj.get_var_traj('C_N'), label='C_N')
            pylab.legend(loc=0)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(0, 1.4)
            pylab.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
            pylab.ylabel('per mRNA (M_P) and\nnuclear PER-TIM complex (C_N)',
                         fontsize='small', horizontalalignment='center')
            for index, tick in enumerate(pylab.arange(0.5, 4.01, 0.5)):
                pylab.text(1.01, float(index)/7.0-0.02, '%.1f'%tick,
                           {'color':'k', 'fontsize': 12}, transform=sp.transAxes)
            pylab.text(1.05, 0.5, 'tim mRNA (M_T)',
                       {'color': 'k', 'fontsize':'small'},
                       horizontalalignment='left',
                       verticalalignment='center',
                       rotation=90, transform=sp.transAxes)
            pylab.title('B')

            sp = pylab.subplot(325)
            pylab.plot(traj.timepoints, traj.get_var_traj('P_t'), label='P_t')
            pylab.plot(traj.timepoints, (traj.get_var_traj('T_t')-1.0)/8.0+1.0,
                       label='T_t')
            pylab.legend(loc=0)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(1.0, 2.0)
            pylab.yticks(pylab.arange(1.0, 2.01, 0.2))
            pylab.ylabel('Total PER protein (P_t)', fontsize='small',
                         horizontalalignment='center')
            for index, tick in enumerate(pylab.arange(1.0, 9.01, 1.0)):
                pylab.text(1.01, float(index)/8.0-0.02, '%.1f'%tick,
                           {'color':'k', 'fontsize': 12}, transform=sp.transAxes)
            pylab.text(1.05, 0.5, 'Total TIM protein (T_t)',
                       {'color': 'k', 'fontsize':'small'},
                       horizontalalignment='left',
                       verticalalignment='center',
                       rotation=90, transform=sp.transAxes)
            pylab.title('C')

        self.assertTrue(net.periodic['stable'],
                        'Failed finding a stable limit cycle.')
        self.assertTrue(net.periodic['tol']<net.periodic['xtol'],
                        'Failed meeting the tolerance criteria.') 
        self.assertAlmostEqual(1.0, net.periodic['period']/period, 4,
                               'Failed for period.')
        for varName, val in net.periodic['stableIC'].items():
            self.assertAlmostEqual(1.0, DDic[varName]/val, 3,
                                   'Failed for %s.'%varName)

    def test_Forced(self):
        """Test that a forced limit cycle can be found"""
        period = 24.0

        net = base_net.copy('test_Forced')
        net.set_var_constant('v_dT', False)
        net.set_var_ic('v_dT', 4.0)
        level = 2
        for t in range(0, 24*level*2, 24):
            net.add_event('Light (%ih)'%t, 'gt(time, %i)'%t, {'v_dT':4.0})
        for t in range(12, 24*level*2, 24):
            net.add_event('Dark (%ih)'%t, 'gt(time, %i)'%t, {'v_dT':2.0})

        for id, val in LDic.items(): net.set_var_ic(id, round(val,2))
        net.set_periodic(period=24.0, xtol=0.00001, level=level, log=True,
                         minVel=0.0001)

        net.Calculate({'M_P':[0, 72]})
        traj = net.trajectory
        
        if plotTraj:
            sp = pylab.subplot(322)
            for var in ['M_P', 'P_t', 'M_T', 'T_t', 'v_dT']:
                pylab.plot(traj.timepoints, traj.get_var_traj(var),
                           label=var)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(0, 10)
            pylab.yticks([0,2,4,6,8,10])
            pylab.ylabel('Total PER and TIM proteins (P_t and T_t)\n'
                         'and their mRNAs (M_P and M_T)', fontsize='small',
                         horizontalalignment='center')
            pylab.legend(loc=0)
            pylab.title('D')

            sp = pylab.subplot(324)
            pylab.plot(traj.timepoints, traj.get_var_traj('M_P'), label='M_P')
            pylab.plot(traj.timepoints, (traj.get_var_traj('M_T')-0.5)/3.5*1.4,
                       label='M_T')
            pylab.plot(traj.timepoints, traj.get_var_traj('C_N'), label='C_N')
            pylab.legend(loc=0)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(0, 1.4)
            pylab.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
            pylab.ylabel('per mRNA (M_P) and\nnuclear PER-TIM complex (C_N)',
                         fontsize='small', horizontalalignment='center')
            for index, tick in enumerate(pylab.arange(0.5, 4.01, 0.5)):
                pylab.text(1.01, float(index)/7.0-0.02, '%.1f'%tick,
                           {'color':'k', 'fontsize': 12}, transform=sp.transAxes)
            pylab.text(1.05, 0.5, 'tim mRNA (M_T)',
                       {'color': 'k', 'fontsize':'small'},
                       horizontalalignment='left',
                       verticalalignment='center',
                       rotation=90, transform=sp.transAxes)
            pylab.title('E')

            sp = pylab.subplot(326)
            pylab.plot(traj.timepoints, traj.get_var_traj('P_t'), label='P_t')
            pylab.plot(traj.timepoints, (traj.get_var_traj('T_t')-1.0)/8.0+1.0,
                       label='T_t')
            pylab.legend(loc=0)
            pylab.xlim(0, 72)
            pylab.xticks([0,12,24,36,48,60,72])
            pylab.xlabel('Time (hr)')
            pylab.ylim(1.0, 2.0)
            pylab.yticks(pylab.arange(1.0, 2.01, 0.2))
            pylab.ylabel('Total PER protein (P_t)', fontsize='small',
                         horizontalalignment='center')
            for index, tick in enumerate(pylab.arange(1.0, 9.01, 1.0)):
                pylab.text(1.01, float(index)/8.0-0.02, '%.1f'%tick,
                           {'color':'k', 'fontsize': 12}, transform=sp.transAxes)
            pylab.text(1.05, 0.5, 'Total TIM protein (T_t)',
                       {'color': 'k', 'fontsize':'small'},
                       horizontalalignment='left',
                       verticalalignment='center',
                       rotation=90, transform=sp.transAxes)
            pylab.title('F')
            
        self.assertTrue(net.periodic['stable'],
                        'Failed finding a stable limit cycle.')
        self.assertTrue(net.periodic['tol']<net.periodic['xtol'],
                        'Failed meeting the tolerance criteria.') 
        self.assertAlmostEqual(1.0, net.periodic['period']/period, 4,
                               'Failed for period.')
        for varName, val in net.periodic['stableIC'].items():
            self.assertAlmostEqual(1.0, LDic[varName]/val, 3,
                                   'Failed for %s.'%varName)

    def test_Nonoscillatory(self):
        """Test that a nonoscillatory network is detected"""
        net = base_net.copy('test_NonOsillatory')

        # Deleting transcription makes the network non-oscillatory
        # approaching the origin (this does not work on unstable fixed pts)
        net.remove_component('M_P transcription')
        net.remove_component('M_T transcription')

        for id, val in LDic.items(): net.set_var_ic(id, round(val,2))

        net.set_periodic(maxfun=250, period=24.0, xtol=0.001, minVel=0.0001,
                         log=True)

        net.Calculate({'M_P':[0, 72]})

        self.assertTrue(net.periodic['stable'],
                        'Failed to approach stable fixed point.')
        self.assertTrue(net.periodic['period']==0.0,
                        'Period (%s) is non-zero, possible oscillations.'%
                        net.periodic['period'])

    def test_deriv_event(self):
        """Test event with derivatives"""
        net = base_net.copy('test_deriv')

        # We're trying to catch an event that would indicate a local maximum
        # of T_2. A known problem is that numerical errors can cause such events
        # to fire near local minima as well.
        net.add_event('T_2_max', 'lt(T_2_deriv_wrt_time,0)', buffer=1e-3)
        traj = Dynamics.integrate(net, [0, 100])

        # Given that events may also fire at local minima, we want to check
        # only that the maxima are found, not necessarily that they are the only
        # events in the trajectory.
        event_times = numpy.array(traj.event_info[0])
        self.assertTrue(numpy.any(abs(event_times - 15.76) < 0.01),
                        'Failed to find maximum at t=15.76, event_times were: '
                        '%s ' % event_times)
        self.assertTrue(numpy.any(abs(event_times - 44.85) < 0.01),
                        'Failed to find maximum at t=44.85, event_times were: '
                        '%s ' % event_times)

    def test_maximum_cost(self):
        """Test cost with maximum"""
        net = base_net.copy('test_max')
        net.add_event('T_2_max', 'lt(T_2_deriv_wrt_time,0)', buffer=1e-3)

        expt = Experiment('Max expt')
        expt.set_fixed_sf({'T_2':1})
        # We need the maxTime here because, in the absence of other data. The
        # trajectory won't get calculated otherwise.
        expt.add_scaled_max('test_max', 'T_2', maxval=5, sigma=0.1, maxTime=50)

        m = Model([expt], [net])
        c = m.cost(m.get_params())

        self.assertAlmostEqual(c, 82.03, 2)

    def test_minimum_cost(self):
        """Test cost with minimum"""
        net = base_net.copy('test_min')
        net.add_event('T_2_min', 'gt(T_2_deriv_wrt_time,0)', buffer=1e-3)

        expt = Experiment('Min expt')
        expt.set_fixed_sf({'T_2':1})
        # This is also a test of restricting range for min/max calculation.
        expt.add_scaled_min('test_min', 'T_2', minval=0.5, sigma=0.1, 
                            minTime=70, maxTime=100)

        m = Model([expt], [net])
        c = m.cost(m.get_params())

        self.assertAlmostEqual(c, 6.00, 2)
        
suite = unittest.makeSuite(test_Periodic)
if __name__ == '__main__':
    unittest.main()
