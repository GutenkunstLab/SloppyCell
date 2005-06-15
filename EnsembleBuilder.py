import Collections
import Residuals
import Model

import sets, scipy, RandomArray, Numeric, copy, pylab, scipy.linalg

class EnsembleBuilder:
    """This was made a class so we don't have to pass around big matrices
    and arrays."""

    def __init__(_,model,currParams):
        _.origParams = copy.deepcopy(currParams)
        _.numParams = len(currParams.keys())
        _.model = copy.deepcopy(model)
        _.usingLogParams = True
	_.keepParamsPositive = False
	_.hessian = scipy.zeros((_.numParams, _.numParams),scipy.Float)
        _.hessFlag = 0
        _.sampMatrix = scipy.zeros((_.numParams, _.numParams),scipy.Float)
        _.matrixV = scipy.zeros((_.numParams, _.numParams),scipy.Float)
        _.singVals = scipy.zeros(_.numParams, scipy.Float)
        _.minSingVal = 0.0 # this is the min sing. val. AFTER we cut off
	_.D = scipy.zeros(_.numParams, scipy.Float)
        _.normSquared = 0.0
        _.cutoff = 1.0e-4  # relative cutoff used in svd
	_.cutoffIndex = 0 # which sing. vals get cut off?
	_.scale = .1 # scale the final step, as it may be still too long
	_.decompFlag = 0
        _.temperature = 1.
        _.ensembleParams = []
        _.ensembleCosts = []
        _.ensembleEigenParams = []
        _.funcAccuracy = 1e-16
        _.DBL_EPSILON = scipy.limits.double_epsilon
        RandomArray.seed(72529486,916423761)


    def BuildEnsemble(_, ensSize):
        currParams = copy.copy(origParams)

	currCost = _.model.Cost(currParams)
        if _.usingLogParams & _.keepParamsPositive :
		print "BuildEnsemble: using log params and params are being kept positive"

	if _.hessFlag == 0:
            if _.usingLogParams :
	    	_.hessian = _.CalcLMHessianInLogSpace()
	    else :
	    	_.CalcHessian(0.001)


        if _.decompFlag == 0:
            _.makeSamplingMatrix()
        _.ensembleParams.append(currParams)
        _.ensembleCosts.append(currCost)
        ensCount = 1
        trialMoves = 0
        while ensCount < ensSize:
            (deltaParams, quadCost) = _.GetTrialMoveAndQuadratic(currParams)
            #print "Trial move is ", deltaParams
	    nextParams = copy.copy(currParams)
            if _.usingLogParams :
	    	nextParams.update(scipy.array(currParams.values())*(scipy.exp(_.scale*scipy.array(deltaParams.values()))))
	    else :
		nextParams.update(scipy.array(currParams.values())+_.scale*scipy.array(deltaParams.values()))

	    tryMove = True
	    if _.keepParamsPositive :
                if min(nextParams[i]) < 0:
                    tryMove = False
	    if tryMove:
	    	newCost = _.model.Cost(nextParams)
		trialMoves += 1
		#print "Quad cost is ", quadCost
		#print "new cost is ", newCost
		if _.AcceptMove(quadCost,newCost-currCost):
			ensCount += 1
			currParams = nextParams
			cost = _.model.Cost(currParams)
			_.ensembleParams.append(currParams)
			_.ensembleCosts.append(cost)
            if trialMoves%200 == 0:
                print "trial moves", trialMoves
                print "ensemble size", ensCount
        print "final trial moves", trialMoves
        return

    def CalcHessianInLogSpace(_) :
	_.hess = _.model.CalcHessianInLogParameters(_.origParams, 
                                           _.funcAccuracy**.25)
        _.hessFlag=1
        return _.hess

    def CalcHessian(_) :
	_.hess = _.model.CalcHessianUsingResiduals(_.origParams, pow(_.funcAccuracy,.25))
        # or could use the following if you want finite differenced Hessian 
        # (be careful
	# because you could step into negative parameter regions)
	# hess = _.model.CalcHessian(currParams, pow(_.funcAccuracy,.25))
	_.hessFlag=1
        return hess

    def CalcLMHessian(_) :
	currParams = _.origParams.__copy__()
        jac, jtj = _.model.GetJandJtJ(currParams)
	return jtj

    def CalcLMHessianInLogSpace(_):
        currParams = _.origParams.__copy__()
        jac, jtj = _.model.GetJandJtJInLogParameters(currParams)

	return jtj

    def makeSamplingMatrix(_, cutoff=scipy.log(2.)):
	## decompose the hessian
        paramsStep = _.origParams.__copy__()

        ## basically need SVD of hessian - singular values and eigenvectors
        ## hessian = u * diag(singVals) * vh
        (u,_.singVals,vh) = scipy.linalg.decomp.svd(_.hessian)

	print "Matrix decomposed. Condition number is ", _.singVals[0]/_.singVals[_.numParams-1]

	## scroll through the singular values and find the ones whose inverses will be
	## huge and set them to zero also, load up the array of singular values that
	## we store
	## cutoff = (1.0/_.singVals[0])*1.0e03
	## double cutoff = _.singVals[0]*1.0e-02
	_.minSingVal = _.cutoff*_.singVals[0]

	foundFlag = False
        for i in range(_.numParams):
            if (_.singVals[i] > _.minSingVal) :
	    	_.D[i] = 1.0/_.singVals[i]
	    else :
		_.D[i] = 1.0/_.minSingVal
		if foundFlag == False :
			_.cutoffIndex = i
			foundFlag = True

	## now fill in the sampling matrix ("square root" of the Hessian)
        ## note that sqrt(D[i]) is taken here whereas Kevin took sqrt(D[j])
        ## this is because vh is the transpose of his PT -JJW
        _.matrixV = vh.copy()
        _.sampMatrix = scipy.transpose(vh)*scipy.sqrt(_.D)
## 	for i in range(_.numParams):
##             for j in range(_.numParams):
##                 _.sampMatrix[i][j] = vh[i][j]*scipy.sqrt(_.D[i])
##                 _.matrixV[i][j] = vh[i][j]

	## Square of the norm of the last row of the sampling matrix
        _.normSquared = Numeric.sum(_.sampMatrix[_.numParams-1]**2)
        _.decompFlag = 1
        return

    def GetTrialMoveAndQuadratic(_,currParams):
        trialMove = currParams.__copy__()
        trialMove.update(scipy.zeros(_.numParams))

	## sample the gaussian in either case

	## draw random numbers r ~ N(0,1)
        randVec = RandomArray.normal(0.,1.,_.numParams)
	# Divide the random vector by an appropriate constant such
	# that the acceptance ratio will always be about exp(-1) if
	# the harmonic approx. is good
	scale = _.cutoffIndex + sum(_.singVals[_.cutoffIndex:])/_.minSingVal
	scale = scipy.sqrt(scale)
	randVec = randVec/scale
	## now rescale by the sampling matrix, including the appropriate
	## reweighting by the control parameter.

	## RIGHT MULTIPLICATION
        trialMove.update(Numeric.dot(_.sampMatrix,randVec))

	## calculate and return the value of the quadratic form
        quadratic = Numeric.dot(Numeric.transpose(trialMove),Numeric.dot(_.hessian,trialMove))

        rmsStep = scipy.sqrt(Numeric.sum(scipy.array(trialMove)**2))/_.numParams
        #print "rmsStep", rmsStep
	return (trialMove, quadratic)

    def AcceptMove(_,quadCost, deltaCost):
        decision = 0
        arg = (-1.*deltaCost)/_.temperature
        if (arg >= 0.0):
            decision = 1
        else:
            p = RandomArray.uniform(0,1)
            if p < scipy.exp(arg):
                decision = 1

        return decision
    
    def PCA(_) :
	# pass in an ensemble (as produced by the ensemble builder)
	# Get out the svd decomposition of the covariance matrix of
	# the log parameters from their mean, which is the principal
	# components analysis
	ensp = _.ensembleParams
	np = len(ensp[0])
	enssize = len(ensp)
	ensarray = scipy.zeros((enssize,np),scipy.Float)
	for i in range(0,enssize) :
		ensarray[i,:] = scipy.array(ensp[i])

	# take logs and subtract off the mean
	logens = scipy.log(ensarray)
	logensmean = scipy.mean(logens,0) # 0 means we take average down columns
	logensmean = scipy.resize(logensmean,(1,np))
	logmeans = scipy.matrixmultiply(scipy.ones((enssize,1),scipy.Float),logensmean)
	displacements = logens - logmeans
	# get covariance matrix
	cov = scipy.matrixmultiply(scipy.transpose(displacements),displacements)
	[u,s,vh] = scipy.linalg.svd(cov)
	return u,s

    def plotHists(_, u, singvals, n) :
	# this is a copy of the same function that was written in Matlab by
	# Ryan. It takes an ensemble, ens, (given back by EnsembleBuilder)
	# projects the displacement of each
	# set of parameters in the ensemble from the mean onto the eigendirections
	# given in u (passed in as columns), and then plots those displacements in
	# a histogram over the range, spread.
	# u could be for example the eigendirections of the Hessian or the PCA
	# directions
	# n is a vector containing which vectors in u you want to use
	# eg. plotHist(ens,p,u,e,[1,3,8],(-3.0,3.0) makes
	# distribution along the u[1][:], u[3][:] and u[8][:] directions
	# Note that unlike other script, we don't sort u or singvals, so for
	# example, u[:,0], singvals[0] may not be the biggest eigenvector/eigenvalue
	# pair
	ensp = _.ensembleParams
	np = len(ensp[0])
	enssize = len(ensp)
	ens = scipy.zeros((enssize,np),scipy.Float)
	for i in range(0,enssize) :
		ens[i,:] = scipy.array(ensp[i])


	logens = scipy.log(ens)
	logensmean = scipy.mean(logens,0) # 0 means we take average down columns

	logensmean = scipy.resize(logensmean,(1,np))
	logmeans = scipy.matrixmultiply(scipy.ones((enssize,1),scipy.Float),logensmean)
	displacements = scipy.log(ens) - logmeans
	pind = 0
	proj = scipy.zeros((enssize,len(n)),scipy.Float)
	for index in n :
		tmp = scipy.matrixmultiply(displacements,u[:,index])
		proj[:,pind] = tmp
		pind = pind+1
	# now for each vector we passed in make a histogram
	bins = 20
	histarrays = []
	for index in range(0,pind) :
		(histarray, bins, patches) = \
		pylab.hist(proj[:,index]*scipy.sqrt(singvals[n[index]]),bins=bins,normed=True)
		histarray = scipy.array(histarray,scipy.Float)
		histarrays.append(histarray)
		# supress the bar histogram output - prefer lines
		pylab.close()
		# print sum(histarray,0)
	for index in range(0,pind) :
		pylab.figure(index+1)
		pylab.plot(bins,histarrays[index],'-')
	pylab.show()

    def GetProjection(_, v, n) :
	# Similar to above except now we just return the ensemble projected onto
	# however many axes we pass in
	ensp = _.ensembleParams
	np = len(ensp[0])
	enssize = len(ensp)
	ens = scipy.zeros((enssize,np),scipy.Float)
	for i in range(0,enssize) :
		ens[i,:] = scipy.array(ensp[i])


	logens = scipy.log(ens)
	logensmean = scipy.mean(logens,0) # 0 means we take average down columns

	logensmean = scipy.resize(logensmean,(1,np))
	logmeans = scipy.matrixmultiply(scipy.ones((enssize,1),scipy.Float),logensmean)
	displacements = scipy.log(ens) - logmeans
	pind = 0
	proj = scipy.zeros((enssize,len(n)),scipy.Float)
	for index in n :
		tmp = scipy.matrixmultiply(displacements,v[:,index])
		proj[:,pind] = tmp
		pind = pind+1
	return proj

    def GetEnsembleTrajs(_,networkcalculation, times) :
	# pass in networkcalculation (an instance of the network)
	# , and a vector of times. 
	# Get out a dictionary indexed by chemical
	# which has the envelope of all the trajectories of all
	# the chemicals over the ensemble. Also passes out the
	# variance at each time point so can compare with the
	# linearized prediction for the variance
	m =  _.model
	clc = networkcalculation

	bestfittraj = {}
	meanchem = {}
	variancechems = {}
	maxtrajchems = {}
	mintrajchems = {}
	varstocalc = {}
	alldvs = clc.dynamicVars.keys() + clc.assignedVars.keys()  # all chemicals
	for var in alldvs :
		maxtrajchems[var] = scipy.zeros((len(times),),scipy.Float)
		mintrajchems[var] = scipy.ones((len(times),),scipy.Float)*1.0e80
		meanchem[var] = scipy.zeros((len(times),),scipy.Float)
		variancechems[var] = scipy.zeros((len(times),),scipy.Float)
		bestfittraj[var] = []
		varstocalc[var] = times


	for var in alldvs :
		# store best fit separately
		for params in _.ensembleParams[0:1] :
			clc.Calculate(varstocalc,params)
			tmptraj = clc.trajectory.getVariableTrajectory(var)
			bestfittraj[var] = scipy.array(tmptraj)

		# this computes trajectories for all parameter sets in the ensemble
		# May want to step through this
		for params in _.ensembleParams :
			clc.Calculate(varstocalc,params)
			tmptraj = clc.trajectory.getVariableTrajectory(var)
			# alltrajs[var].append(tmptraj)
			meanchem[var] = meanchem[var] + scipy.array(tmptraj)
			variancechems[var] = variancechems[var] + scipy.array(tmptraj)**2 # an elementwise square
			for i in range(0,len(times)) :
				maxtrajchems[var][i] = max(maxtrajchems[var][i],tmptraj[i])
				mintrajchems[var][i] = min(mintrajchems[var][i],tmptraj[i])

	for kys in meanchem.keys() :
		meanchem[kys] = meanchem[kys]/len(_.ensembleParams)
		variancechems[kys] = variancechems[kys]/len(_.ensembleParams) - meanchem[kys]**2

	return bestfittraj,meanchem,variancechems,maxtrajchems,mintrajchems
