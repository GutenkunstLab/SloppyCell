import scipy

# see VdmPairwise for definition of r2PC12 and dtMatPC12
# 1. make cluster instance
# PC12cluster = cluster.HierarchicalClustering(range(48), lambda x,y: getDist(x,y,r2PC12/dtMatPC12))
# 2. set desired linake method (see cluster.py for different options)
# PC12cluster.setLinkageMethod('single')
# 3. getlevel returns the clusters of items with requested distance of one another.
#    You must run this before you can ask for the topology
# PC12cluster.getlevel(0.1)
# 4. get topology (formatted as a nested tuple
# PC12topo = PC12cluster.topo()
# 5. convert topo into a permutation list
# PC12permList = FlattenIt(PC12topo,r2PC12/dtMatPC12)
# 6. get permutation matrix from permutation list
# PC12permMat = makePermutationMatrix(permList)
# 7. permute your jacobian
# PC12permJac = scipy.dot(jacPC12,PC12permMat)

def getDist(ind1,ind2,distMat):
    """
    distMat is a distance matrix.  distMat[i,j] == distance from i to j
    """
    return distMat[ind1,ind2]

def FlattenIt(x, distMat):
    """
    Takes topo output from hierarchical clustering and orders the list at the bottom
    level to correspond to the optimal permutation.
    """
    if isinstance(x, int):
	return [x]
    a = FlattenIt(x[0], distMat)
    b = FlattenIt(x[1], distMat)
#    print a, b
    d00 = distMat[a[0],b[0]]
    d01 = distMat[a[0],b[-1]]
    d10 = distMat[a[-1],b[0]]
    d11 = distMat[a[-1],b[-1]]
    dMin = min(d00, d01, d10, d11)
    if d00 == dMin:
	a.reverse()
    elif d01 == dMin:
	a.reverse()
	b.reverse()
    elif d11 == dMin:
        b.reverse()
    a.extend(b)
#    print a
    return a

def makePermutationMatrix(permList):
    """
    Takes a list defining the permutation and makes the appropriate matrix.
    """
    permList = scipy.array(permList)
    n = len(permList)
    if 0 not in permList:
        permList = permList - 1
    permMat = scipy.zeros((n,n),'d')
    for ii, jj in enumerate(permList):
        permMat[ii,jj] = 1.
    return scipy.transpose(permMat)


def normColumns(mat1):
    """
    Normalizes each column of a matrix.  So far just for display purposes.
    """
    n = len(mat1[0])
    mat2 = mat1.copy()
    for ii in range(n):
        mat2[:,ii] = mat2[:,ii]/scipy.sqrt(scipy.dot(mat2[:,ii],mat2[:,ii]))
    return mat2

def colorAlignColumns(mat1):
    """
    Multiplies columns by +/- 1 to make them have positive dot product with neighboring
    columns.  So far just for display purposes.
    """
    n = len(mat1[0])
    mat2 = mat1.copy()
    for ii in range(n-1):
        if scipy.dot(mat2[:,ii],mat2[:,ii+1]) < 0.:
            mat2[:,ii+1] = -1.*mat2[:,ii+1]
    return mat2
