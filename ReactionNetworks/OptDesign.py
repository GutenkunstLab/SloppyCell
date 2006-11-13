import scipy, copy
import SloppyCell.Utility
load = SloppyCell.Utility.load
save = SloppyCell.Utility.save
import SloppyCell.ReactionNetworks.Dynamics as Dynamics

try:
    import SloppyCell.Plotting as Plotting
except ImportError:
    pass

def setup(paramfile,calcobject,senstrajfile,jtjfile) :
    """ Set up the quantities necessary to run the optimal design
    algorithms. NOTE: This function needs to be called first 
    before any of the optimal design functions can be called.
    
    paramfile: the name of a pickled file containing the 
    best fit parameters in KeyedList format
    calcobject: the calculation object for which we are doing the 
    optimal design. (Note that in general, may be searching a 
    design over many different calculations, but here we only 
    consider one. Thus, we set design_sentraj equal to senstraj)
    senstrajfile: the name of the file containing the pickled 
    sensitivity trajectory for the calculation, calcobject, 
    for the set of parameters in paramfile.
    jtjfile: the name of the file containing the pickled Fisher
    Information Matrix (J^t J) for the current set of data and
    for the parameters in paramfile. 
    NOTE: The derivatives computed for J^tJ need to be with respect
    to the *log* of the parameters
    """
    import OptDesign as v
    v.curp = load(paramfile)
    v.jtj = load(jtjfile)
    v.clc = calcobject 
    v.senstraj = load(senstrajfile)
    v.design_senstraj = v.senstraj
    v.p_names_ordered = v.curp.keys()
    v.jtjdict = {}
    for pindex1,pname1 in enumerate(v.p_names_ordered) :
        for pindex2,pname2 in enumerate(v.p_names_ordered) :
            v.jtjdict[(pname1,pname2)] = v.jtj[pindex1][pindex2]

    v.ovvarnames = v.clc.optimizableVars.keys() 
    v.jtjtrunc = scipy.zeros((len(v.ovvarnames),len(v.ovvarnames)),scipy.float_)

    # The number of optimizable variables for the calculation we are 
    # considering might be less than the number of parameters for the 
    # whole model. We are only working with this calculation so we 
    # need to trim down the J^t J (Fisher information) matrix 
    # accordingly
    for pindex1,pname1 in enumerate(v.ovvarnames) :
        for pindex2,pname2 in enumerate(v.ovvarnames) :
            v.jtjtrunc[pindex1][pindex2] = v.jtjdict[(pname1,pname2)]


def make_sens_traj(calcobject,params,times,senstrajfilename):
    """ Make the sensitivity trajectory for the calculation
    calcoject (same as in setup(...) above). 
    params: parameters as a KeyedList, sensitivity traj is 
    calculated at these parameters (should be same as in paramfile 
    in setup(...) above)
    times: the timepoints in the sensitivity trajectory (1-d array)
    senstrajfilename: the file to save the sensitivity trajectory to

    Note that if times is very finely spaced, the 
    sensitivity trajectory will need a lot of storage space """
    senstraj = Dynamics.integrate_sensitivity(calcobject, times, params, 1.0e-6)
    save(senstraj,senstrajfilename)

def design_over_chems(chemnames,designchemnames,logprior=1.0e20) :
    """ 
    chemnames = list of unmeasurable chemicals
    designchemnames = list of measurable chemicals
    logprior = prior on params, e.g. logprior = log(1000.0) means
    parameter standard deviation will be less than a factor of 1000.0
    
    Out of the list chemnames, find the best chemical and 
    best time point, that most reduces the integrated variance
    over designchemnames """
    times = design_senstraj.timepoints    
    trunc_times = [times[i] for i in scipy.arange(0,len(times),1)]    
    best_change = 0.0 # the change should always be negative    
    best_chem = "None"
    best_time = "None"
    for dchemname in designchemnames :
        print "On design chemical ", dchemname    
        for t in trunc_times :
            sensvect_design = get_sens_vect(dchemname,t)
            # NOTE: assuming a 10% error on the measurement --- use 10% of the 
            # maximum value in the trajectory
            maxval = max(design_senstraj.get_var_traj(dchemname)) + 1.0
            sensvect_design = sensvect_design/(.1*maxval)
            intvar_change = integrated_var_change(chemnames,sensvect_design,logprior)
            tot_change = 0.0    
            for id in chemnames :
                tot_change = tot_change + intvar_change[id]
            if tot_change < best_change :
                best_change = tot_change
                best_chem = dchemname
                best_time = t
    return best_change, best_chem, best_time

def design_over_single_variance(sensvect,designchemnames,logprior=1.0e20) :
    """
    sensvect : a sensitivity vector (length = # of params) of 
    unmeasurable quantity of interest
    designchemnames : list of measurable chemicals
    
    sensvect could be the sensitivity of a single chemical at a
    single timepoint; then can use method get_sens_vect (see elsewhere
    in this file) to compute this sensitivity vector. In that
    case we are designing over the species variance at that single point
    """
    times = senstraj.timepoints    
    trunc_times = [times[i] for i in scipy.arange(0,len(times),5)]    
    best_change = 0.0 # the change should always be negative    
    best_chem = "None"
    best_time = "None"
    for dchemname in designchemnames :
        for t in trunc_times :
            sensvect_design = get_sens_vect(dchemname,t)
            var_change = single_variance_change(sensvect,sensvect_design,logprior)
            if var_change < best_change :
                best_change = var_change
                best_chem = dchemname
                best_time = t
    return best_change, best_chem, best_time

def variances(chemnames,logprior=1.0e20) :
    """ chemnames : list of chemical names for which the 
    variance at all timepoints will be computed
    logprior : prior on parameters. logprior = log(1000.0)
    means params allowed to vary by about a factor of 1000.0
    return values : 
    times: times of the trajectory
    bestfit: a dictionary of best fit trajectories (keys are entries in chemnames)
    var: a dictionary of variances (keys are entries in chemnames)
    """
    #senstraj = load('EndogenousEGFR3T3sensNoPriors')
    times = senstraj.timepoints    
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(
        len(jtjtrunc),len(jtjtrunc)))    
    var = {}
    bestfit = {}    
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    for name in chemnames :
        var[name] = []    
        bestfit[name] = []    
        chemindex = senstraj.key_column.get(name)    
        index1sens = senstraj.key_column.get((name,first))
        index2sens = senstraj.key_column.get((name,last))
        sensarray_this_chem = copy.copy(senstraj.values[:,index1sens:(index2sens+1)])
        # Turn sensitivities into sensitivities with respect to log parameters
        for j, pname in enumerate(ovvarnames) :
            sensarray_this_chem[:,j] = sensarray_this_chem[:,j]*curp.get(pname)
        
        tmp = scipy.dot(sensarray_this_chem,jtjinv)
        for i in range(len(tmp[:,0])) :
            var[name].append(scipy.dot(tmp[i,:],sensarray_this_chem[i,:]))
            
        bestfit[name] = senstraj.values[:,chemindex]
        var[name] = scipy.asarray(var[name])    
    return times, bestfit, var

def variances_log_chems(chemnames,logprior=1.0e20) :
    """ Same as above except the variances are now on the 
    logs of the chemicals trajectories.
    """
    #senstraj = load('EndogenousEGFR3T3sensNoPriors')
    times = senstraj.timepoints    
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(
        len(jtjtrunc),len(jtjtrunc)))    
    var = {}
    bestfit = {}    
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    for name in chemnames :
        var[name] = []    
        bestfit[name] = []    
        chemindex = senstraj.key_column.get(name)    
        index1sens = senstraj.key_column.get((name,first))
        index2sens = senstraj.key_column.get((name,last))
        sensarray_this_chem = copy.copy(senstraj.values[:,index1sens:(index2sens+1)])
        traj_this_chem = copy.copy(senstraj.values[:,chemindex])    
        for j, pname in enumerate(ovvarnames) :
            sensarray_this_chem[:,j] = sensarray_this_chem[:,j]*curp.get(pname)
        # need to scale each row by 1/chemvalue to mimic a derivative w.r.t. 
        # log chemicals. Add a small value to chemvalue to avoid divide by zero
        for i in range(len(times)) :
            sensarray_this_chem[i,:] = sensarray_this_chem[i,:]/(traj_this_chem[i]+1.0e-6)

        tmp = scipy.dot(sensarray_this_chem,jtjinv)
        for i in range(len(tmp[:,0])) :
            var[name].append(scipy.dot(tmp[i,:],sensarray_this_chem[i,:]))
            
        bestfit[name] = senstraj.values[:,chemindex]
        var[name] = scipy.asarray(var[name])    
    return times,bestfit,var

def single_variance(sensvect,logprior=1.0e20) :
    """ Get the variance for a single function of parameters 
    that has a sensitivity vector sensvect. Useful for looking at
    variances in parameter combinations, or simple functions of
    parameters. Note that if we are concerned with ratios and 
    products of parameters, it's often best to consider sensvect
    as a sensitivity w.r.t. log parameters """
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(
        len(jtjtrunc),len(jtjtrunc)))    
    tmp = scipy.dot(jtjinv,sensvect)
    var = scipy.dot(sensvect,tmp)
    return var

def variance_change(chemnames,sensvect_design,logprior=1.0e20) :
    """
    chemnames : list of chemical names at which we will look 
    at variance
    sensvect_design : the sensitivity vector (one by no. params array) at
    the new design point.
    returns : (times, varchange) 
    the times and the change in variances at those times (should
    be negative) for each of the chemicals in chemnames, after the 
    addition of the new timepoint. varchange is a dictionary 
    indexed by entries in chemnames.
    """
    times = senstraj.timepoints    
    n = len(jtjtrunc)    
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(n,n))

    #sensvect_design = scipy.resize(sensvect_design,(n,1))
    jtjinv_design = scipy.dot(jtjinv,sensvect_design)
    #jtjinv_design = scipy.resize(jtjinv_design,(n,1))     # want a column vector    

    denominator = 1.0 + scipy.dot(sensvect_design,jtjinv_design)
    
    varchange = {}
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    for name in chemnames :
        varchange[name] = []    
        chemindex = senstraj.key_column.get(name)    
        index1sens = senstraj.key_column.get((name,first))
        index2sens = senstraj.key_column.get((name,last))
        sensarray_this_chem = copy.copy(senstraj.values[:,index1sens:(index2sens+1)])
        for j, pname in enumerate(ovvarnames) :
            sensarray_this_chem[:,j] = sensarray_this_chem[:,j]*curp.get(pname)
        
        product = scipy.dot(sensarray_this_chem,jtjinv_design)
        # this product is a number of timepoints by one vector, we need to 
        # square each element for the final formula
        varchange[name] = -scipy.asarray(product**2/denominator)

    return times, varchange

def single_variance_change(sensvect,sensvect_design,logprior=1.0e20) :
    """
    sensvect : given a single function f(p) of parameters, this is the
    derivative w.r.t. each of the parameters (in log parameters). For
    ratios or products of rate constants, f(p) is a linear function
    sensvect_design : the sensitivity vector of the new point in the 
    design you wish to add
    returns: the variance change of the quantity f(p), given the 
    addition of the new data point, with sensitivity vector sensvect_design.
    """
    n = len(jtjtrunc)    
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(n,n))

    jtjinv_design = scipy.dot(jtjinv,sensvect_design)
    denominator = 1.0 + scipy.dot(sensvect_design,jtjinv_design)
    product = scipy.dot(sensvect,jtjinv_design)
    return -product**2/denominator

def get_sens_vect(chemname,time) :
    """ get a sensitivity vector for a chemical "chemname" at a
    time, time """
    tindex = design_senstraj._get_time_index(time,1.0e-4)
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    index1sens = design_senstraj.key_column.get((chemname,first))
    index2sens = design_senstraj.key_column.get((chemname,last))
    sens_vect = copy.copy(
            design_senstraj.values[tindex,index1sens:(index2sens+1)])
    for j, pname in enumerate(ovvarnames) :
        sens_vect[j] = sens_vect[j]*curp.get(pname)
    return sens_vect

def get_sens_array(chemname) : 
    """ get an array of sens_vects for all the times the chemical is defined
    and convert to log sensitivities """    
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    chemindex = design_senstraj.key_column.get(chemname)    
    index1sens = design_senstraj.key_column.get((chemname,first))
    index2sens = design_senstraj.key_column.get((chemname,last))
    sensarray_this_chem = copy.copy(
            design_senstraj.values[:,index1sens:(index2sens+1)])
    for j, pname in enumerate(ovvarnames) :
        sensarray_this_chem[:,j] = sensarray_this_chem[:,j]*curp.get(pname)
    return sensarray_this_chem

def integrated_var_change(chemnames,sensvect_design,logprior=1.0e20) :
    times, varchange = variance_change(chemnames,sensvect_design,logprior)
    int_varchange = {}
    for name in varchange.keys() :
        int_varchange[name] = scipy.integrate.simps(varchange[name],times)
    
    return int_varchange

def var_change_weighted(weights,chemnames,sensarray_design,logprior=1.0e20) :
    """ This is similar to var_change except now we pass in a sensarray
    instead of sensvect --- this is a matrix of sensvects aligned rowwise.
    Row i will be multiplied by sqrt(weights[i]) where sum(weights)=1 and 
    each weight is a number between zero and one. We will return the 
    change in variance for all the chemicals in chemnames """
    # we use the formula (Sherman-Woodbury-Morrison)
    # (A+UV^t)^(-1) = A^(-1) - A^(-1)*U*(I + V^T*A^(-1)*U)^(-1)*V^t*A^(-1)
    # where U = V and V^t = W^(1/2)*sensarray_design  
    times = senstraj.timepoints    
    ntimes = len(times)    
    k,n = sensarray_design.shape    
    jtjinv = scipy.linalg.inv(jtjtrunc+1.0/logprior**2*scipy.eye(n,n))
    Vt = scipy.zeros((k,n),scipy.float_)
    for i in range(k) :
        Vt[i,:] = scipy.sqrt(weights[i])*sensarray_design[i,:]
    design_jtjinv = scipy.dot(Vt,jtjinv)
    #jtjinv_design = scipy.resize(jtjinv_design,(n,1))     # want a column vector    

    denominator = scipy.eye(k,k) + \
            scipy.dot(design_jtjinv,scipy.transpose(Vt))
    inv_denom = scipy.linalg.inv(denominator)

    varchange = {}
    optvarkeys = clc.optimizableVars.keys()
    first = optvarkeys[0]
    last = optvarkeys[-1]
    for name in chemnames :
        varchange[name] = []    
        chemindex = senstraj.key_column.get(name)    
        index1sens = senstraj.key_column.get((name,first))
        index2sens = senstraj.key_column.get((name,last))
        sensarray_this_chem = copy.copy(senstraj.values[:,index1sens:(index2sens+1)])
        for j, pname in enumerate(ovvarnames) :
            sensarray_this_chem[:,j] = sensarray_this_chem[:,j]*curp.get(pname)
        
        product = scipy.dot(design_jtjinv,
                scipy.transpose(sensarray_this_chem))
        # each column vector of this matrix has to be dotted through the
        # denominator matrix --- each column is a different time point
        for j in range(ntimes) :    
            quadprod = scipy.dot(product[:,j],inv_denom)    
            quadprod = scipy.dot(quadprod,product[:,j])    
            varchange[name].append(-quadprod)
        varchange[name] = scipy.asarray(varchange[name])
    
    return times, varchange

def integrated_var_change_weighted(weights,chemnames,sensarray_design,logprior=1.0e20) :
    times, varchange = var_change_weighted(weights,chemnames,sensarray_design,
            logprior)
    intvarchange = {}
    for name in varchange.keys() :
        intvarchange[name] = scipy.integrate.simps(varchange[name],times)
    return intvarchange

def weight_cost(weights,chemnames,sensarray_design,logprior=1.0e20) :
    """ For this cost function we're going to assume unconstrained 
    variables are being passed in, so we need to convert them to 
    a range between 0 and 1. The sum of the weights should also = 1 """
    weights0to1 = weights_trans(weights)
    # now weights lie between 0 and 1    
    weights0to1 = weights0to1/scipy.sum(weights0to1) # this makes sure 
    # weights sum up to 1.
    intvarchange = integrated_var_change_weighted(weights0to1,chemnames,
            sensarray_design,logprior)
    cost = 0.0    
    for n in intvarchange.keys() :
        cost = cost + intvarchange[n]
    return cost

def weights_trans(weights) :
    wtrans = (scipy.sin(weights)+1.0)/2.0
    return wtrans

def weights_inv_trans(transweights) :
    w = scipy.arcsin(2.0*transweights-1.0)
    return w

def minimize_weight_cost(weights,chemnames,sensarray_design,logprior=1.0e20) :
    """
    weights : a vector of positive numbers with length the same as the number of 
    rows of sensarray_design. The weights should sum to 1 
    chemnames: a list of unmeasurable chemical names over which we wish
    to design experiments 
    sensarray_design: an array of sensitivities of measurable chemicals
    or just an array of sensitivity vectors, each row a different 
    sensitivity vector 
    logprior : prior on parameters. logprior = log(1000.0) allows parameters
    to fluctuate by a factor of 1000 """      
    weights_trans = scipy.arcsin(2.0*weights-1.0)    
    # maxiter may need to be increased if convergence is not apparent
    # or if the number of weights is increased
    w = scipy.optimize.fmin(weight_cost,weights_trans,maxiter = 10000,
            args=(chemnames,sensarray_design,logprior))
    woptnotnormed = (scipy.sin(w)+1.0)/2.0
    wopt = woptnotnormed/scipy.sum(woptnotnormed)
    return woptnotnormed,wopt

def plot_variances(chemnames,logprior,scale=1.0,return_var = False) :
    """ 
    chemnames: list of chemical names 
    logprior: prior on params. logprior = log(1000.0) means parameters
    allowed to fluctuate by a factor of 1000 """
    times, bestfit, var = variances(chemnames,logprior)
    for key in bestfit.keys() :
        Plotting.figure()    
        Plotting.plot(times,bestfit[key]/scale)
        Plotting.hold(True)    
        Plotting.plot(times,bestfit[key]/scale + scipy.sqrt(var[key])/scale,'r--')
        Plotting.plot(times,bestfit[key]/scale - scipy.sqrt(var[key])/scale,'r--')
        Plotting.title(key,fontsize=16)
        Plotting.xlabel('time (minutes)',fontsize=16)
        Plotting.ylabel('number of molecules',fontsize=16)
        xtics = Plotting.gca().get_xticklabels()
        ytics = Plotting.gca().get_yticklabels()
        Plotting.setp(xtics,size=16)
        Plotting.setp(ytics,size=16)
        #Plotting.axis([0.0,40.0,-.01,1.2e4])    
    Plotting.show()
    if return_var :
        return times, bestfit, var

def plot_variances_log_chems(chemnames,logprior) :
    """
    chemnames: list of chemical names
    logprior: prior on params
    Plots the standard deviation of the chemicals when the variance
    is computed using logs of the chemical trajectories. This 
    makes sure the final plots do not have best_fit+-stddev that
    do not become negative """
    times, bestfit, var = variances_log_chems(chemnames,logprior)
    for key in bestfit.keys() :
        Plotting.figure()    
        Plotting.plot(times,bestfit[key])
        Plotting.hold(True)    
        Plotting.plot(times,bestfit[key]*scipy.exp(scipy.sqrt(var[key])),'r-')
        Plotting.plot(times,bestfit[key]*scipy.exp(-scipy.sqrt(var[key])),'r-')
        Plotting.title(key,fontsize=14)
        Plotting.xlabel('time')
        Plotting.ylabel('arb. units')
        #Plotting.axis([0.0,40.0,-.01,1.2e4])    
    Plotting.show()

def plot_variance_newpoint(chemnames,sensvect_design,logprior=1.0e20,
        return_data = True) :
    """
    chemnames: list of chemical names
    sensvect_design: a sensivity vector of a quantity that is 
    measurable
    This will plot the old and new variances of the chemicals in 
    chemnames, given a new measurement that has sensitivity vector
    sensvect_design
    """
    times,bestfit,var = variances(chemnames,logprior)
    times,varchange = variance_change(chemnames,sensvect_design,logprior)    

    for key in bestfit.keys() :
        Plotting.figure()    
        Plotting.plot(times,bestfit[key])
        Plotting.hold(True)    
        Plotting.plot(times,bestfit[key] + scipy.sqrt(var[key]),'r-')
        Plotting.plot(times,bestfit[key] - scipy.sqrt(var[key]),'r-')
    
        Plotting.plot(times,bestfit[key] + scipy.sqrt(var[key]+varchange[key]),'k--')
        Plotting.plot(times,bestfit[key] - scipy.sqrt(var[key]+varchange[key]),'k--')
        
        Plotting.title(key,fontsize=14)
        Plotting.xlabel('time')
        Plotting.ylabel('arb. units')
        Plotting.axis([0.0,40.0,-.01,1.2e4])    
    Plotting.show()
    if return_data :
        newvar = {}
        for ky in var.keys() :
            newvar[ky] = var[key] + varchange[key]
        return times,bestfit,newvar

def plot_variance_newweights(weights,chemnames,sensarray_design,logprior=1.0e20,scale=1.0,return_data = True) :
    """
    weights : a proposed set of weights for each of the row vectors in 
    sensarray_design
    chemnames : a list of chemicals for which we will plot the variance
    logprior : as before 
    This will plot the old and new variances on chemnames, similar to
    above.
    NOTE: the weights that are passed in do not necessarily have to sum to
    one. e.g. if the weights are normalized such that max(weights) = 1, then
    by scaling all the weights by 1/sigma, you are then assuming that
    the most accurate measurement has an error of size sigma. sigma for 
    example could be 20% of the maximum value of a trajectory.
    """

    times,bestfit,var = variances(chemnames,logprior)
    times,varchange = var_change_weighted(weights,chemnames,sensarray_design,logprior)    

    for key in bestfit.keys() :
        Plotting.figure()    
        Plotting.plot(times,scale*bestfit[key])
        Plotting.hold(True)    
        Plotting.plot(times,scale*bestfit[key] + scale*scipy.sqrt(var[key]),'r-')
        Plotting.plot(times,scale*bestfit[key] - scale*scipy.sqrt(var[key]),'r-')
    
        Plotting.plot(times,scale*bestfit[key] + scale*scipy.sqrt(var[key]+varchange[key]),'k--')
        Plotting.plot(times,scale*bestfit[key] - scale*scipy.sqrt(var[key]+varchange[key]),'k--')
        
        Plotting.title(key,fontsize=14)
        Plotting.xlabel('time')
        Plotting.ylabel('arb. units')
        Plotting.axis([0.0,40.0,-.01,1.2e4])    
    Plotting.show()

    if return_data :
        newvar = {}
        for ky in var.keys() :
            newvar[ky] = var[key] + varchange[key]
        return times,bestfit,newvar

def plot_variances_subplot(chemnames,logprior) :
    times, bestfit, var = variances(chemnames,logprior)
    nallplots = len(chemnames)
    # 9 at a time    
    nfigs = nallplots/9 # integer division -- no fractional part
    for figno in range(1,nfigs+1) :
        Plotting.figure()
        for i in range(0,9) :    
            Plotting.subplot(3,3,i+1)    
            chemind = i+(figno-1)*9    
            Plotting.plot(times,bestfit[chemnames[chemind]])
            Plotting.hold(True)
            Plotting.plot(times,bestfit[chemnames[chemind]] 
                    + scipy.sqrt(var[chemnames[chemind]]),'r-')
            Plotting.plot(times,bestfit[chemnames[chemind]] 
                    - scipy.sqrt(var[chemnames[chemind]]),'r-')
            
            yt = Plotting.yticks()
            Plotting.axis([0,100.0,yt[0],yt[-1]])    
            Plotting.title(chemnames[chemind])
            Plotting.xlabel('time')
            Plotting.ylabel('arb. units')
            xt = Plotting.xticks()
            Plotting.xticks([xt[0],xt[-1]])
        
        Plotting.savefig('./figs/variance_wt_'+i.__str__()+'.ps')    
    Plotting.show()


#def fix_sf():
    # make sure scale factors get computed --- easiest way is
    # to compute the cost
#    print "cost is ", m.cost(curp)
#    sfs = m.internalVars['scaleFactors']
#    for exptname in sfs.keys() :
#        fixeddict = sfs[exptname]    
#        m.exptColl[exptname].set_fixed_sf(fixeddict)
    # just check
#    print "cost is now", m.cost(curp)

def reduce_size(array,skipsize) :
    """ reduce_size takes an array of dimension m,n and 
    returns an array with every skipsize row sampled.
    """
    size = array.shape
    newsize = len(scipy.arange(0,size[0],skipsize))    
    if len(size) == 1 : # a vector 
        newvect = scipy.zeros((newsize,),scipy.float_)
        for iind,i in enumerate(scipy.arange(0,size[0],skipsize)) :
            newvect[iind] = array[i]
        return newvect
    elif len(size) == 2 : # an array
        newarray = scipy.zeros((newsize,size[1]),scipy.float_)
        for iind,i in enumerate(scipy.arange(0,size[0],skipsize)) :
            newarray[iind] = array[i]
        return newarray
