from __future__ import division
from builtins import range
from past.utils import old_div
from . import OptimizeSumDets as OSD
try:
    import cluster
except ImportError:
    pass
from . import clusterScripts
from . import VdmPairwise as Vdm
import scipy

def getPermList(hess):
    Dmat = Vdm.calcDmat(hess)
    lambdaMat, r2Mat =  Vdm.getLambdaR2(hess,Dmat)
    ctMat, dtMat = Vdm.calcR2mat_Terms(hess, lambdaMat)
    clHess = cluster.HierarchicalClustering(list(range(len(hess))), lambda x,y: clusterScripts.getDist(x,y,old_div(r2Mat,dtMat)))
    clHess.getlevel(0.1)
    permList = clusterScripts.FlattenIt(clHess.topo(),old_div(r2Mat,dtMat))
    return permList
