import logging
logger = logging.getLogger('script')
import sys
import time
import traceback

import scipy

from SloppyCell.ReactionNetworks import *
Network.full_speed()
Dynamics.global_rtol = 1e-6
try:
    import SloppyCell.ReactionNetworks.RunInParallel as Par
    Par.statement_to_all_workers('Network.full_speed()')
    Par.statement_to_all_workers('Dynamics.global_rtol = %g' 
                                 % Dynamics.global_rtol)
except (ImportError, NameError):
    logger.warn('Failure importing RunInParallel.')

import Experiments as Expts
import Calculations as Calcs
from min5_wFinal import min5_wFinal_params

import Nets
reload(Nets)
m = Model_mod.Model([Expts.ErkMekTraverse2EGF.expt,
           Expts.Raf1LandrethEGF.expt,
           Expts.Erk1HYaoEGF.expt,
           Expts.ErkMekTraverse2NGF.expt,
           Expts.RasGreen1NGF.expt,
           Expts.Rap1YorkNGF.expt,
           Expts.Erk1HYaoNGF.expt,
           Expts.RasGreen1EGF.expt,
           Expts.Raf1BRafLandrethNGF.expt,
           Expts.ErkMekTraverseHERover.expt,
           ],
          [Nets.EGFstim100,
           Nets.NGFstim50,
           Nets.EGFstim30,
           Nets.NGFstim100,
           Nets.EGFRx50_EGFstim100,
           ])

params = KeyedList(min5_wFinal_params)
params.setOrder(m.get_params().keys())

# This is a prior that runs (with 95% probability) from value/prior_range to
#  value*prior_range
#prior_range = 1e3
#for key, value in params.items():
#    res = Residuals.PriorInLog('%s_prior' % key, key, scipy.log(value),
#                               scipy.log(scipy.sqrt(prior_range)))
#    m.AddResidual(res)

# We add events to ensure that we catch it when species concentrations fall
#  below their corresponding Michaelis Menten constants.
#for net in m.get_calcs().values():
#    for rxn in net.reactions:
#        if isinstance(rxn, Reactions.MichaelisMentenReaction):
#            for substrate, stoich in rxn.stoichiometry.items():
#                if stoich == -1:
#                    break
#            for Km in rxn.parameters:
#                if Km.startswith('Km'):
#                    break
#            event_name = '%s_%s_declining' % (substrate, Km)
#            if not net.events.has_key(event_name):
#                net.add_event(event_name, 'lt(%s, %s)'% (substrate, Km))

# Starting not from zero, but from tiny amount, helps the integrator get
#  started much more cleanly at tighter tolerances.
#for net in m.get_calcs().values():
#    for id in net.dynamicVars.keys():
#        if net.get_var_ic(id) == 0:
#            net.set_var_ic(id, 1e-6)

m.compile()

import inspect, types
def make_bumping_func(in_func):
    if inspect.ismethod(in_func):
        # If in_func is a method of an object, we need to capture the object
        #  it is a method of and ensure that the object gets passed in at the
        #  'self' position.
        self = in_func.im_self
        def func(*args, **kwargs):
            return in_func.im_func(self, *args, **kwargs)
        func_name = in_func.im_func.func_name
    else:
        func = in_func
        func_name = in_func.func_name

    def bumping_func(*args, **kwargs):
        rtol_storage = rtol = Dynamics.global_rtol
        succeeded = False
        X = None

        while (not succeeded) and (rtol >= 1e-12):
            try:
                ret = func(*args, **kwargs)
                succeeded = True
            except Utility.SloppyCellException, X:
                rtol /= 10
                Dynamics.global_rtol = rtol
                Par.statement_to_all_workers('Dynamics.global_rtol = %g' % rtol)

        if X:
            # Make sure we restore rtol to where it started.
            Dynamics.global_rtol = rtol_storage
            Par.statement_to_all_workers('Dynamics.global_rtol = %g'
                                         % rtol_storage)

        if succeeded:
            if X:
                logger.debug('Difficult evaluation of wrapped function "%s".' 
                             % func_name)
                logger.debug('Succeeded with rtol = %g.' % rtol)
                logger.debug('Arguments were: %s, %s.'
                             % (str(args), str(kwargs)))
            return ret
        else:
            logger.debug('Difficult evaluation of wrapped function "%s".'
                         % func_name)
            logger.debug('Never succeeded with rtol up to %g.' % rtol)
            raise X

    if inspect.ismethod(func):
        types.MethodType(bumping_func, obj, obj.__class__)

    return bumping_func

costs, times = [], []
for ii in range(4):
    start = time.time()
    costs.append(m.cost(params))
    times.append(time.time()-start)
print 'Costs: %s' % ', '.join(['%7.4f'%c for c in costs]) 
print 'Took:  %s' % ', '.join(['%7.2f'%t for t in times])

m_log = m.copy()
print 'Uniform over B free energy: %.2f' % m.free_energy(params, T=1)
for expt_id, expt in m_log.get_expts().items():
    for group in expt.get_sf_groups():
        curr_sf = m.internalVars['scaleFactors'][expt_id][group[0]]
        expt.set_sf_prior(group, 'gaussian in log sf',
                          (scipy.log(curr_sf), scipy.log(scipy.sqrt(1e3))))

print
costs, times = [], []
for ii in range(4):
    start = time.time()
    costs.append(m_log.cost(params))
    times.append(time.time()-start)
print 'Costs: %s' % ', '.join(['%7.4f'%c for c in costs]) 
print 'Took:  %s' % ', '.join(['%7.2f'%t for t in times])
print 'Gaussian over logB free energy: %.2f' % m_log.free_energy(params, T=1)

#e = Utility.load('final_culled_1e3_prior_ensemble.bp')
#jtjs = []
#for p in e[::10]:
#    j = m.jacobian_log_params_sens(scipy.log(p))
#    jtj = scipy.dot(scipy.transpose(j), j)
#    jtjs.append(jtj)
#Utility.save(jtjs, 'more_jtjs.bp')

#gold = m.free_energy(params, 1)
#jtj = Utility.load('jtj.best_energy.bp')
#u_orig,v = Utility.eig(jtj)
#jtj3 = Utility.load('PCA_hessian.j.bp')
#u3,v3 = Utility.eig(jtj3)
#u = u_orig
##u = u3
#from scipy import *
#jtj_p = sum([u[ii]*outer(v[:,ii], v[:,ii]) for ii in range(48)], axis=0)
#samp_mat = Ensembles._sampling_matrix(jtj3)
#acceptances = []
#diffs = []
#for ii in range(500):
#    move = Ensembles._trial_move(samp_mat)
#    pnew = scipy.asarray(params)*scipy.exp(move)
#    diff = Ensembles._quadratic_cost(move, jtj)
#    gnew = m.free_energy(pnew, 1)
#    diff = gnew-gold
#    diffs.append(diff)
#    accepted = Ensembles._accept_move(diff, 1)
#    acceptances.append(accepted)

#ej, gj, rj = Utility.load('ens.log_gaussian.jtj.6.tight.bp')
#elem = scipy.random.randint(len(ej))
#print 'Starting from log_gaussian.jtj.6.tight element %i' % elem
#params = ej[elem]
#jtj = Utility.load('mean_inv_jtj.j.bp')
#
#m.free_energy = make_bumping_func(m.free_energy)
#m.GetJandJtJInLogParameters = make_bumping_func(m.GetJandJtJInLogParameters)
#print 'Current free energy: %.2f' % m.free_energy(params, T=1)
#try:
#    Dynamics.logger.setLevel(Dynamics.logging.CRITICAL)
#    Dynamics.daskr.logger.setLevel(Dynamics.daskr.logging.CRITICAL)
#    save_to = 'ens.temp.bp'
#    Ensembles.logger.setLevel(Ensembles.logging.DEBUG)
#    e, c, r = Ensembles.ensemble_log_params(m, scipy.array(params), 
#                                            hess = jtj,
#                                            steps = 1000,
#                                            step_scale=1.5,
#                                            sing_val_cutoff=1e-3,
#                                            recalc_hess_alg = False,
#                                            skip_elems = 100,
#                                            save_to = save_to,
#                                            save_hours = 1.)
#except Exception, X:
#    print 'Uncaught exception!'
#    tb = traceback.format_exception(sys.exc_type, sys.exc_value, 
#                                    sys.exc_traceback)
#    print ''.join(tb)
#    print 'Current pararmeters are %s.' % str(m.get_params())
#
#skip = int(1.*len(ej)/50.)
#jtjs = []
#for p in ej[::skip]:
#    start = time.time()
#    j, jtj = m.GetJandJtJInLogParameters(scipy.log(p))
#    jtjs.append(jtj)
#    taken = time.time() - start
#    print 'Took %7.2f' % taken
#    Utility.save(jtjs, 'ens.log_gauss.jtj.6.jtjs.bp')

#times = []
#m.attach_observer('print', Observers.CostPrinter())
#for ii in range(1):
#    start = time.time()
#    j, jtj = m.GetJandJtJInLogParameters(scipy.log(params))
#    times.append(time.time()-start)
#print 'JtJs'
#print 'Took:  %s' % ', '.join(['%7.2f'%t for t in times])

#m.cost = make_bumping_func(m.cost)

#ens, costs, r = Utility.load('ens.tight.1.bp')
#diff_c = scipy.diff(costs)
#attempts_with_change = scipy.arange(1, len(costs))[diff_c != 0]
#jtjs = [JtJ]
#successes = [0]
#attempts = [0]
##for ii in scipy.linspace(0, len(ens)-1, 50):
#for success in 2**scipy.arange(0, 19):
#    success = int(success)
#    attempt = attempts_with_change[success-1]
#    params = ens[attempt]
#    j, jtj = m.GetJandJtJInLogParameters(scipy.log(params))
#    print 'Done with jtj number %i, %i.' % (success, attempt)
#    jtjs.append(jtj)
#    successes.append(success)
#    attempts.append(attempt)
#Utility.save((successes, attempts, jtjs), 'jtjs.ens.tight.1.2.bp')

#for var in net.optimizableVars.keys():
#    net.set_var_optimizable(var, False)
#net.set_var_optimizable('krbEGF', True)

#for net in m.get_calcs().values():
#    net.compile()
#    net.compile_c()
#
#costs, times = [], []
#for ii in range(4):
#    start = time.time()
#    costs.append(m.cost(params))
#    times.append(time.time()-start)
#print 'Costs: %s' % ', '.join(['%7.4f'%c for c in costs]) 
#print 'Took:  %s' % ', '.join(['%7.2f'%t for t in times])

#old_net = net.copy('old')
#
## def res_func(time, dynamicVars, yprime, ires, constants, residual):
#import scipy.weave
#mod = scipy.weave.ext_tools.ext_module('Brown_res_ext')
#time = 0.0
#dynamicVars = scipy.ones(32, scipy.float_)
#yprime = scipy.ones(32, scipy.float_)
#constants = scipy.ones(49, scipy.float_)
#residual = scipy.ones(32, scipy.float_)
#ires = 0
#res = scipy.weave.ext_tools.ext_function('res_func', net.res_func_c_body,
#                                         ['time', 'dynamicVars', 'yprime', 
#                                          'ires', 'constants', 'residual'])
##res.customize.add_support_code(res_code)
#mod.add_function(res)
#mod.compile()
#import Brown_res_ext
#net.res_function = Brown_res_ext.res_func

#for ii in range(int(1e4)):
#    old_net.res_func(time, dynamicVars,yprime, ires, constants,residual)
#net.res_func(time, dynamicVars,yprime, ires, constants,residual)

#Utility.disable_warnings()
#m.cost(hard_params)
#Utility.enable_warnings()
#m.cost(hard_params)

#J, JtJ = m.GetJandJtJInLogParameters(scipy.log(params))
#jac = scipy.array(J.values())
#permList = Vandermonde.decompHess.getPermList(JtJ)
#permMat = Vandermonde.clusterScripts.makePermutationMatrix(permList)
#Plotting.pcolor(Vandermonde.Vdm.colorAlignColumns(Vandermonde.Vdm.normColumns(scipy.dot(jac,permMat))))
#Plotting.pcolor(Vandermonde.Vdm.colorAlignColumns(Vandermonde.Vdm.normColumns(scipy.dot(scipy.transpose(permMat),scipy.dot(JtJ,permMat)))))

#m.cost = lambda p: Ensembles._quadratic_cost(scipy.log(p) - scipy.log(p_array),
#                                             jtj)

#jtj_new = m.GetJandJtJInLogParameters(scipy.log(params))[1]
#Utility.save(jtj, 'jtj.no_prior.bp')

#ens, costs, r = Utility.load('finally.1.bp')
#jtjs = []
#for ii in scipy.linspace(0, len(ens)-1, 100):
#    ii = int(ii)
#    params = ens[ii]
#    jtj = None #m.GetJandJtJInLogParameters(scipy.log(params))[1]
#    print 'Done with jtj number %i.' % ii
#    jtjs.append(jtj)
#
#Utility.save(jtjs, 'finally.jtjs.bp')

#ens, costs, ratio = Utility.load('ens.tight.1.bp')
#cost_jtj = sample_jtj = Utility.load('jtj.tight.bp')
#
#from scipy import log
#
#successful_steps = [0] + scipy.compress(scipy.diff(costs) != 0, range(1, len(ens)+1)).tolist()
#to_try = [0]+[int(res) for res in scipy.logspace(0, 4, 100)]
#to_try = scipy.sort(scipy.unique(to_try))
#ratios = []
#steps_taken_l = []
#steps_tried_l = []
#for steps_taken in to_try:
#    steps_tried = successful_steps[steps_taken]
#    params = ens[steps_tried]
#    if steps_tried != 0:
#        cost_jtj = m.GetJandJtJInLogParameters(scipy.log(params))[1]
#    m.cost = lambda p: Ensembles._quadratic_cost(log(p)-log(params), cost_jtj) 
#    e, c, ratio = Ensembles.ensemble_log_params(m, scipy.array(params), 
#                                                hess = sample_jtj, 
#                                                steps=int(1e2))
#    print steps_tried, steps_taken, ratio
#    steps_tried_l.append(steps_tried)
#    steps_taken_l.append(steps_taken)
#    ratios.append(ratio)
#    Utility.save((steps_tried_l, steps_taken_l, ratios), 'using_quads.bp')

