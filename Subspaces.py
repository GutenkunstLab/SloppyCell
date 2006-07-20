import scipy

def subspace_angles(A, B):
    """
    Return the principle angles and vectors between two subspaces. 

    See Bjorck and Golub, "Numerical Methods for Computing Between Linear 
           Subspaces" _Mathematics_of_Computation_, 27(123), 1973, pp. 579-594
    Or, for a more understandable exposition, Golub and Van Loan, 
           _Matrix_Computations_ (3rd ed.) pp. 603-604
    """
    A, B = scipy.asarray(A), scipy.asarray(B)
    if A.shape[0] != B.shape[0]:
        raise ValueError, 'Input subspaces must live in the same dimensional '\
                'space.'

    # Get orthogonal bases for our two subspaces.
    QA, QB = scipy.linalg.orth(A), scipy.linalg.orth(B)

    M = scipy.matrixmultiply(scipy.transpose(scipy.conjugate(QA)), QB)
    Y, C, Zh = scipy.linalg.svd(M)

    U = scipy.matrixmultiply(QA, Y)
    V = scipy.matrixmultiply(QB, scipy.transpose(scipy.conjugate(Zh)))
    return scipy.arccos(C), U, V

def test():
    # This is the example from Bjorck and Golub (1973)
    m, p = 26, 13

    A = scipy.zeros((m, p), scipy.Float)
    for col in range(p):
        A[col*(m/p):(col+1)*(m/p), col:col+1] = scipy.ones((m/p, 1))
    A *= 1/scipy.sqrt(m/p)

    B = scipy.zeros((m, p), scipy.Complex)
    for col in range(0, p):
        B[:, col:col+1] = scipy.transpose(scipy.mat((-1 + 2j*scipy.arange(m)/(m+1))**col))

    theta, U, V = subspace_angles(A, B)
    assert scipy.round(scipy.cos(theta), 4) == \
            scipy.array([ 1.    ,  0.9982,  0.9981,  0.9903,  0.9899,  0.9765,
                         0.9628,  0.9415, 0.9176,  0.8701,  0.7637,  0.0608,  
                         0.0156]) 
    assert scipy.round(scipy.sin(theta), 4) == \
            scipy.array([ 0.    ,  0.0594,  0.0609,  0.1388,  0.1418,  0.2157,
                         0.2701,  0.337 , 0.3975,  0.4928,  0.6456,  0.9982,  
                         0.9999])
    print 'Comparison with Bjorck and Golub (1973) example passed.'

    # Example from p 605 of Matrix Computations
    A = [[1,2],[3,4],[5,6]]
    B = [[1,5],[3,7],[5,-1]]
    theta, U, V = subspace_angles(A, B)
    assert scipy.round(scipy.cos(scipy.real(theta)), 3) ==\
            scipy.array([1.000, 0.856])
    
def plot_test():
    A = scipy.rand(3, 2) - 0.5
    B = scipy.rand(3, 1) - 0.5
    theta, U, V = subspace_angles(A, B)
    
    import os, tempfile
    file_names = []
    for vects in [A, B, U, V]:
        fd, name = tempfile.mkstemp()
        file_names.append(name)
        f = os.fdopen(fd, 'w')
        lines = []
        for column in scipy.transpose(vects):
            column = scipy.asarray(column)
            normalized = column/scipy.sqrt(scipy.dot(column, column))
            lines.append(' '.join([str(0) for elem in column]))
            lines.append(' '.join([str(elem) for elem in normalized]))
            lines.append('')
            lines.append('')
        f.write(os.linesep.join(lines))
        f.close()
    
    gnuplot_command = 'splot "%s" lw 10 t "A", "%s" lw 10 t "B", '\
            '"%s" lw 8 t "U", "%s" lw 4 t "V"' % tuple(file_names)

    print theta * 180/scipy.pi
    print gnuplot_command
