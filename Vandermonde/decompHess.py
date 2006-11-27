import OptimizeSumDets as OSD
try:
    import cluster
except ImportError:
    pass
import clusterScripts
import VdmPairwise as Vdm
import scipy

def getPermList(hess):
    Dmat = Vdm.calcDmat(hess)
    lambdaMat, r2Mat =  Vdm.getLambdaR2(hess,Dmat)
    ctMat, dtMat = Vdm.calcR2mat_Terms(hess, lambdaMat)
    clHess = cluster.HierarchicalClustering(range(len(hess)), lambda x,y: clusterScripts.getDist(x,y,r2Mat/dtMat))
    clHess.getlevel(0.1)
    permList = clusterScripts.FlattenIt(clHess.topo(),r2Mat/dtMat)
    return permList
