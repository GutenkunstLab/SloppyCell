import copy
import os
import sets

import scipy
try:
    import scipy.misc.limits as limits
except ImportError:
    limits = scipy.limits

import SloppyCell._lsodar as _lsodar
import SloppyCell.Utility as Utility

_msgs = {2: "Integration successful.",
         3: "Integration successful. Root found.",
         -1: "Excess work done on this call (perhaps wrong Dfun type).",
         -2: "Excess accuracy requested (tolerances too small).",
         -3: "Illegal input detected (internal error).",
         -4: "Repeated error test failures (internal error).",
         -5: "Repeated convergence failures (perhaps bad Jacobian or tolerances).",
         -6: "Error weight became zero during problem.",
         -7: "Internal workspace insufficient to finish (internal error)."
         }

_rwork_vars = {'hu': 10,
              'tcur': 12,
              'tolsf': 13,
              'tsw': 14,}
_iwork_vars = {'nst': 10,
              'nfe': 11,
              'nje': 12,
              'nqu': 13,
              'imxer': 15,
              'lenrw': 16,
              'leniw': 17,
              'mused': 18,
              }

class odeintrException(Utility.SloppyCellException):
    pass

def odeintr(func, y0, t, args=(), Dfun=None, full_output=0, ml=0, mu=0, rtol=1e-6, atol=1e-12, tcrit=None, h0=0.0, hmax=0.0, hmin=0.0, ixpr=0, mxstep=500, mxhnil=0, mxordn=12, mxords=5, printmessg=0, root_term = [], root_func=None, int_pts=False, insert_events=False, return_derivs = None):

    """Integrate a system of ordinary differential equations.

    Description:

      Solve a system of ordinary differential equations Using lsodar from the
      FORTRAN library odepack.
      
      Solves the initial value problem for stiff or non-stiff systems
      of first order ode-s:
           dy/dt = func(y,t0,...) where y can be a vector.

      Additionally, will search for roots of the optional vector function 
      root_func. If root_term[index of root found] = 1, the integration will
      terminate at that point.

      Note that if roots are sought, the return sequence is yout, tout, i_root,
      t_root, y_root. yout and tout will contain all points at which roots were
      found. i_root contains the indices of the roots found, t_root contains 
      when they were found, and y_root contains the y-values when they were 
      found.

    Inputs:

      func -- func(y,t0,...) computes the derivative of y at t0.
      y0   -- initial condition on y (can be a vector).
      t    -- a sequence of time points for which to solve for y.  The intial 
              value point should be the first element of this sequence.
      args -- extra arguments to pass to function.
      Dfun -- the gradient (Jacobian) of func (same input signature as func).
      full_output -- non-zero to return a dictionary of optional outputs as
                     the second output.
      printmessg -- print the convergence message.
      root_term -- a sequence declaring whether or not to terminate the
                   integration when a given root is found. If root_term[i] is 
                   True the integration will terminate upon finding root i.
      root_func -- root_func(y, t0) computes the function whose roots are sought
                   Note: If roots are sought, the return values are yout, tout,
                   t_root, y_root, i_root. t_root, y_root, i_root are the time
                   values, y values, and indices of the roots found. The times
                   at which roots were found are inserted into yout and tout.
                   

    Outputs: (y, {tout, t_root, y_root, i_root, infodict,})

      y -- a rank-2 array containing the value of y in each row for each
           desired time in t (with the initial value y0 in the first row).

      infodict -- a dictionary of optional outputs:
        'hu'    : a vector of step sizes successfully used for each time step.
        'tcur'  : a vector with the value of t reached for each time step.
                  (will always be at least as large as the input times).
        'tolsf' : a vector of tolerance scale factors, greater than 1.0,
                  computed when a request for too much accuracy was detected.
        'tsw'   : the value of t at the time of the last method switch
                  (given for each time step).
        'nst'   : the cumulative number of time steps.
        'nfe'   : the cumulative number of function evaluations for each
                  time step.
        'nje'   : the cumulative number of jacobian evaluations for each
                  time step.
        'nrfe'  : the cumulative number of root function evaluations for
                  each timestep
        'nqu'   : a vector of method orders for each successful step.
        'imxer' : index of the component of largest magnitude in the
                   weighted local error vector (e / ewt) on an error return.
        'lenrw' : the length of the double work array required.
        'leniw' : the length of integer work array required.
        'mused' : a vector of method indicators for each successful time step:
                  1 -- adams (nonstiff)
                  2 -- bdf (stiff)

    Additional Inputs:

      *** Not tested yet!
      ml, mu -- If either of these are not-None or non-negative, then the
                Jacobian is assumed to be banded.  These give the number of
                lower and upper non-zero diagonals in this banded matrix.
                For the banded case, Dfun should return a matrix whose
                columns contain the non-zero bands (starting with the
                lowest diagonal).  Thus, the return matrix from Dfun should
                have shape len(y0) x (ml + mu + 1) when ml >=0 or mu >=0
      ***          
      rtol -- The input parameters rtol and atol determine the error
      atol    control performed by the solver.  The solver will control the
              vector, e, of estimated local errors in y, according to an
              inequality of the form
                   max-norm of (e / ewt) <= 1
              where ewt is a vector of positive error weights computed as
                   ewt = rtol * abs(y) + atol
              atol can be either a vector the same length as y or a scalar.
              rtol must be a scalar.
      tcrit -- a vector of critical points (e.g. singularities) where
               integration care should be taken.

       (For the next inputs a zero default means the solver determines it).

      h0 -- the step size to be attempted on the first step.
      hmax -- the maximum absolute step size allowed.
      hmin -- the minimum absolute step size allowed.
      ixpr -- non-zero to generate extra printing at method switches.
      mxstep -- maximum number of (internally defined) steps allowed
                for each integration point in t.
      mxhnil -- maximum number of messages printed.
      mxordn -- maximum order to be allowed for the nonstiff (Adams) method.
      mxords -- maximum order to be allowed for the stiff (BDF) method.

      int_pts -- If True, every step the integrator took is output. The return
                sequence is then yout, tout, [info_dict]

      return_derivs -- If True, the return sequence becomes:
                          yout, tout, [ie, ye, te,] dout, [info_dict].
                          dout is an array of the 1st time derivatives of
                          the variables.
    """

    # Take care of our passed in functions. We need to wrap them in lambda's to 
    #  pass on the extra arguments
    usefunc = lambda t, y: func(y, t, *args)
    if len(root_term) > 0 and root_func is not None:
        useg = lambda t, y: root_func(y, t, *args)
    else:
        useg = lambda x: x
    if Dfun is not None:
        usejac = lambda t, y: Dfun(y, t, *args)
        Dfun_def = True
    else:
        usejac = lambda x: x
        Dfun_def = False

    neq = len(y0)
    ng = len(root_term)

    # Now we take care of the tolerance settings
    if not scipy.isscalar(rtol):
        raise ValueError, 'rtol must be a scalar.'

    if not scipy.isscalar(atol) and len(atol) != neq:
        raise ValueError, 'atol must be a scalar or a vector with same length as number of equations.'

    # itol is 2 if atol is a vector, 1 otherwise.
    itol = 2 - scipy.isscalar(atol)

    itask = 1
    # Do we have critical points?
    if tcrit:
        itask = 4
    if int_pts:
        # By default, the integrator won't stop at tout if it's in intermediate
        #  output mode. To make sure we return results for the requested times,
        #  we add them to our critical times list.
        if tcrit:
            tcrit_set = sets.Set(tcrit)
            tcrit_set.union_update(t[1:])
            tcrit = list(tcrit_set)
            tcrit.sort()
        else:
            tcrit = t[1:]
        itask = 5

    # Is our jacobian banded?
    banded_jac = (ml > 0) or (mu > 0)

    # Set jacobian type
    jt = {(True, False): 1,
          (False, False): 2,
          (True, True): 4,
          (False, True): 5}[Dfun_def, banded_jac]


    # Make the work arrays
    lmat = max(neq**2 + 2, (2*ml + mu + 1)*neq + 2)
    lrn = 20 + neq*(mxordn+4) + 3*ng
    lrs = 20 + neq*(mxords+4) + lmat + 3*ng 
    lrw = max(lrn, lrs)
    rwork = scipy.zeros(lrw, scipy.float_)

    liw = 20 + neq
    iwork = scipy.zeros(liw, scipy.int_)

    rwork[4] = h0
    rwork[5] = hmax
    rwork[6] = hmin
    iwork[0] = ml
    iwork[1] = mu
    iwork[4] = ixpr
    iwork[5] = mxstep
    iwork[6] = mxhnil
    iwork[7] = mxordn
    iwork[8] = mxords

    istate = 1 # First call
    t0, tindex = t[0], 1
    tout, yout = [t0], [y0]
    t_root, y_root, i_root = [], [], []

    # Handle derivative requests
    if return_derivs:
        dout = [copy.copy(usefunc(t0, y0))]

    info_dict = dict([(key, []) for key in
                      (_rwork_vars.keys() + _iwork_vars.keys())])
    tcrit_ii = 0
    while tindex < len(t):
        twanted = t[tindex]
        if (itask == 4 or itask == 5) and tcrit_ii < len(tcrit):
            rwork[0] = tcrit[tcrit_ii]
            twanted = min(tcrit[tcrit_ii], twanted)

        y, treached, istate, jroot = \
                _lsodar.dlsodar(usefunc, copy.copy(y0), t0, twanted, 
                                itol, rtol, atol, 
                                itask, istate, rwork, iwork, 
                                usejac, jt, 
                                useg, ng)

        if istate < 0:
            # Problem!
            print _msgs[istate]
            print "Run with full_output = 1 to get quantitative information."
            outputs = (scipy.array(yout), tout, t_root, y_root, i_root)
            if full_output:
                outputs = outputs + (info_dict,)
            raise odeintrException(_msgs[istate],outputs)
        else:
            if printmessg:
                print _msgs[istate]

            # If we need to record this point.
            if treached == t[tindex] or itask == 5 or\
               (istate == 3 and insert_events):
                yout.append(y)
                tout.append(treached)
                if full_output:
                    for key, index in _rwork_vars.items():
                        info_dict[key].append(rwork[index])
                    for key, index in _iwork_vars.items():
                        info_dict[key].append(iwork[index])
                    info_dict['message'] = _msgs[istate]

                # Collect derivatives
                if return_derivs:
                    t_ask = treached * (1- limits.double_resolution)
                    dky, iflag = _lsodar.dintdy(t_ask, 1, rwork[20+ng:],
                                                neq)
                    dout.append(dky)

            # If we reached our goal, move on to the next point.
            if treached == t[tindex]:
                tindex += 1
                    
            # If we reached a critical point, move to the next
            if (itask == 4 or itask == 5) and treached == tcrit[tcrit_ii]:
                tcrit_ii += 1
                # If we're out of critical points, drop back to normal 
                #  integration mode
                if tcrit_ii == len(tcrit):
                    itask = 1

            # We've found a root
            if istate == 3:
                jroot = list(jroot)
                outputs = (scipy.array(yout), tout, t_root, y_root, i_root)
                if full_output:
                    outputs = outputs + (info_dict,)
                if jroot.count(1) > 1:
                    print 'Multiple roots found at a single point!?! '\
                            'jroot is %s' % jroot
                    raise odeintrException(_msgs[istate],outputs)
                elif jroot.count(1) == 0:
                    print 'LSODAR claimed root found, but jroot is empty. '\
                            'jroot is %s' % jroot
                    raise odeintrException(_msgs[istate],outputs)

                # Which root did we hit?
                crossed = jroot.index(1)

                t_root.append(treached)
                y_root.append(copy.copy(y))
                i_root.append(crossed)

                if root_term[crossed] == 1:
                    break

    # Process outputs
    outputs = (scipy.array(yout), tout, t_root, y_root, i_root)
    if return_derivs:
        outputs = outputs + (scipy.array(dout),)

    if full_output:
        outputs = outputs + (info_dict,)

    return outputs
