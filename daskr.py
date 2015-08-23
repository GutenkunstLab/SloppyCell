import logging
logger = logging.getLogger('daskr')

import scipy

import SloppyCell._daskr as _daskr
import SloppyCell.Utility as Utility

# we won't need the psol function, but it's a required argument for DDASKR,
# so I'm defining a dummy function here. if the user does not pass a jac
# function, we still need to pass a dummy to daskr
def dummy_func():
    pass

# These messages correspond to the values of IDID on lines 992-1069 of ddaskr.f
_msgs = {0: "Unexepected status code. (IDID=0)",
         1: "Integration successful. Call the code again to continue the \
integration another step in the direction of TOUT. (IDID=1)",
         2: "Integration successful. Define a new TOUT and call the code \
again.  TOUT must be different from T.  You cannot change the \
direction of integration without restarting. (IDID=2)",
         3: "Integration successful. Define a new TOUT and call the code again.\
TOUT must be different from T.  You cannot change the direction of \
integration without restarting. (IDID=3)",
         4: "Reset INFO(11) = 0 and call the code again to begin the \
integration. (If you leave INFO(11) > 0 and INFO(14) = 1, you may \
generate an infinite loop.) In this situation, the next call to \
DDASKR is considered to be the first call for the problem in that all \
initializations are done. (IDID=4)",
         5: "Call the code again to continue the integration in the direction \
of TOUT. You may change the functions Ri defined by RT after a return \
with IDID = 5, but the number of constraint functions NRT must remain \
the same. If you wish to change the functions in RES or in RT, then \
you must restart the code. (IDID=5)",
         -1: "The code has taken more steps than max_steps (default = \
500). If you want to continue, call the function again with a higher \
value for max_steps. (IDID=-1)",
         -2: "The error tolerances RTOL, ATOL have been increased to values \
the code estimates appropriate for continuing. You may want to change \
them yourself. If you are sure you want to continue with relaexed \
error tolerances, set INFO(1) = 1 and call the code again. (IDID=-2)",
         -3: "A solution component is zero and you set the corresponding \
component of ATOL to zero.  If you are sure you want to continue, you \
must first alter the error criterion to use positive values of ATOL \
for those components corresponding to zero solution components, then \
set INFO(1) = 1 and call the code again. (IDID=-3)",
         -5: "Your JAC routine failed with the Krylov method. Check for errors \
in JAC and restart the integration. (IDID=-5)",
         -6: "Repeated error test failures occurred on the last attempted step \
in DDASKR.  A singularity in thesolution may be present.  If you are \
absolutely certain you want to continue, you should restart the \
integration.  (Provide initial values of Y and YPRIME which are \
consistent.) (IDID=-6)",
         -7: "Repeated convergence test failures occured on the last attempted \
step in DASKR. An inaccurate or ill-conditioned Jacobian or \
preconditioner may be the problem. If you are absolutely certain you \
want to continue, you should restart the integration. (IDID=-7)",
         -8: "The matrix of partial derivatives is singular, with the use of \
the direct methods. Some of your equations may be redundant. DDASKR \
cannot solve the problem as stated. It is possible that the redundant \
equations could be removed, and then DDASKR could solve the problem. \
It is also possible that a solution to your problem either does not \
exist or is not unique. (IDID=-8)",
         -9: "DDASKR has multiple convergence test failures, preceded by \
multiple error test failures, on the last attempted step. It is \
possible that your problem is ill-posed and cannot be solved using \
this code. Or, there may be a discontinuity or a singularity in the \
solution. If you are absolutely certain you want to continue, you \
should restart the integration. (IDID=-9)",
         -10: "DDASKR had multiple convergence test failures because IRES was \
equal to -1. If you are absolutely certain you want to continue, you \
should restart the integration. (IDID=-10)",
         -11: "There was an unrecoverable error (IRES=-2) from RES inside the \
nonlinear system solver. Determine the cause before trying again. \
(IDID=-11)",
         -12: "DDASKR failed to compute the initial Y and YPRIME vectors. \
This could happed because the initial approximation to Y or YPRIME \
was not very good, or because no consistent values of these vectors \
exist. The problem could also be caused by an inaccurate or singular \
iteration matrix, or a poor preconditioner. (IDID=-12)",
         -13: "There was an unrecoverable error encountered inside your PSOL \
routine. Determine the cause before trying again. (IDID=-13)",
         -14: "The Krylove linear solver failed to achieve convergence. \
(IDID=-14)",
         -33: "You cannot continue the solution of this problem. An attempt to \
do so will result in your run being terminated. (IDID=-33)",
         }

redir = Utility.Redirector()

class daeintException(Utility.SloppyCellException):
    pass


def daeint(res, t, y0, yp0, rtol, atol, nrt = 0, rt = None, jac = None, 
            tstop=None, intermediate_output=False, ineq_constr=False,
            calculate_ic=False, var_types=None, redir_output=True,
            max_steps=500, rpar=None, hmax=None, max_timepoints=scipy.inf):
    """Integrate a system of ordinary differential algebraic equations with or
    without events.

    Description:

      Solve a system of ordinary differential algebraic equations using the
      FORTRAN library DDASKR.
      
      Uses the backward differentiation formulas of orders one through five
      to solve a system of the form G(T,Y,Y') = 0 where Y can be a vector.

      Values for Y and Y' at the initial time must be given as input.  These
      values should be consistent, that is, if t0, y0, yp0 are the given
      initial values, they should satisfy G(t0,y0,yp0) = 0.  However, if
      consistent values are not known, in many cases you can have DDASKR
      solve for them if the appropriate option in the info array is set.

      Normally, DDASKR solves the system from some intial t0 to a specified
      tout.  This wrapper takes a list of times t as input, and then solves
      the system for the range of times.  In the interval mode of operation,
      only the solution values at time points in t will be returned.
      Intermediate results can also be obtained easily by specifying info(3).

      While integrating the given DAE system, DDASKR also searches for roots
      of the given constraint functions Ri(T,Y,Y') given by RT. If DDASKR
      detects a sign change in any Ri(T,Y,Y'), it will return the intermediate
      value of T and Y for which Ri(T,Y,Y') = 0.

      Caution: If some Ri has a root at or very near the initial time, DDASKR
      may fail to find it, or may find extraneous roots there, because it does
      not yet have a sufficient history of the solution.
      
    Inputs:

      res  -- res(t, y, ...) computes the residual function for the
              differential/algebraic system to be solved.
      t    -- a sequence of time points for which to solve for y.  The intial 
              value of the dependent variable(s) y should correspond to the
              first time in this sequence.
      y0   -- this array must contain the initial values of the NEQ solution
              components at the initial time point (note NEQ is defined as the
              length of this array).
      yp0  -- this array must contain the initial values of the NEQ first
              derivatives of the solution components at the initial time point.
      jac  -- this is the name of a routine that you may supply that relates to
              the Jacobian matrix of the nonlinear system that the code must
              solve at each time step.
      rtol -- array of relative error tolerances (must be non-negative).
      atol -- array of absolute error tolerances (must be non-negative).
      rt   -- this is the name of the subroutine for defining the vector
              Ri(T, Y, Y') of constraint functions.
      nrt  -- the number of constraint functions Ri(T, Y, Y').
      jac  -- the name of a routine that you may supply (optionally) that
              relactes to the Jacobian matrix of the nonlinear system that the
              code must solve at each time step. If no Jacobian is passed then
              the system will automatically use a numerical Jacobian.
      args -- parameter for passing extra information to the wrapper.
      tstop -- sometimes it is not possible to integrate past some point tstop
              because the equation changes there or it is not defined past
              tstop.  if tstop is set the integrator will not integrate past
              tstop.  NOTE:  If you select a tstop, the integrator may skip
              some points in your list of time points, t.  This is because it
              may pass tstop before hitting all the points in the list of times
              for an integration that is progressing quickly.
      intermediate_output -- True to have all output returned.
      ineq_constr -- True to have inequality constraint checking on.
              If True, the code will check for non-negativity in Y during the
              integration (not during initial condition calculations).
      calculate_ic -- True to have the DDASKR automatically calculate
             consistent initial conditions.  Given Y_d, the code will
             calculate Y_a and Y'_d (note: you must also pass var_types
             to indicate which variables are algebraic and which are
             differential).  We do not use the functionality of daskr
             that allows calculation of Y given Y'.
      var_types -- this list tells whether each variable in the system is
              differential or algebraic. Vor a variable Yi, var_types[i] = +1 if
              Yi is differential, and var_types[i] = -1 if it is algebraic.
      redir_output -- False to have print statements from DASKR printed to the
              standard output.  If redir_output is True, the output will be
              redirected.
      max_steps -- set this if you want the integrator to be able to take more
              than 500 steps between time points when you're not in
              "intermediate output" mode.  If intermediate_output = True then
              every time step is considered an output point so there are never
              multiple steps between time points.  If intermediate_output =
              False (default) then the integrator will stop after taking 500
              steps unless max_steps is set > 500.  The default value for
              max_steps is 500.
      rpar -- If not None, this sequence is passed to the function being
              evaluated and can be used to pass additional arguments.
      hmax -- If not None, this sets an upper limit to the step sizes that 
              ddaskr will take during the integration. Setting this can help
              integration for difficult systems, but it can also slow things
              down.
         
    Outputs:  (yout, tout, ypout, t_root, y_root, i_root)

      yout -- the numerically calculated list of solution points corresponding
              to the times in tout.
      tout -- the list of time points corresponding to yout.
      ypout -- the numerically calculated list of solution points corresponding
      to the times in tout.
      t_root -- the time of a terminating event.
      y_root -- the system state at time t_root.
      i_root -- tells which event(s) fired and how they fired.
                If i_root(i) = 0, then event i did not fire. If nonzero,
                i_root(i) shows the direction of the sign change in Ri.
                i_root(i) = 1  means Ri changed from negative to positive.
                i_root(i) = -1 means Ri changed from positive to negative.
      
    """

    if scipy.isscalar(rpar):
	raise ValueError('rpar must be a sequence.')

    if rpar is None or len(rpar) == 0:
        # rpar needs to have at least length 1 or daskr will screw up arguments.
        rpar = scipy.empty(1)

    y = y0
    yp = yp0
    ires = 0

    # We calculate the number of equations in the residual function from the
    # length of the dependent variable array.
    # The number of event functions are passed to daeint as nrt.
    neq = len(y0)
    # ipar is used to hold neq and the length of rpar.
    ipar = [neq, len(rpar)]

    # The size of the work arrays normally depends on the options we select.
    # Since info(11)=0 always and info(8)=0 always in this wrapper, we
    # can always set BASEr, BASEi, and MAXORD as follows:
    MAXORD = 5
    BASEr = 60+max(MAXORD+4,7)*neq+3*nrt
    BASEi = 40
    # Since info(5) = 0, add neq*2 to LRW
    LRW = BASEr + (neq**2)
    LIW = BASEi + 2*neq

    # Make the work arrays
    # Making these int32 and float64 avoids crashes on 64-bit.
    rwork = scipy.zeros(LRW, scipy.float64)
    iwork = scipy.zeros(LIW, scipy.int32)

    # The info array is used to give the daskr code more details about how we
    # want the problem solved.
    # Set all the options for info, as well as related options.
    info = scipy.zeros(20, scipy.int32)

    # Now we go through the options one by one and set them according to the
    # inputs. Note that some options are not used and kept at '0'.

    # Is this the first call for the problem? (yes)
    info[0] = 0
    
    # Are both error tolerances rtol and atol scalars? (no, we use arrays of
    # tolerances)
    info[1] = 1

    # Now we take care of the tolerance settings
    if scipy.isscalar(rtol) or len(rtol) != neq:
        raise ValueError('rtol must be an array or a vector with same length '
                         'as number of equations.')
    if scipy.isscalar(atol) or len(atol) != neq:
        raise ValueError('atol must be an array or a vector with same length '
                         'as number of equations.')

    # Do you want the solution at intermediate steps (and not just at the
    # specified points)?
    if intermediate_output == False:
        info[2] = 0
    elif intermediate_output == True:
        info[2] = 1
    else:
        raise ValueError('int_output must be True or False.')

    # Can the integration be carried out without any restrictions on the
    # independent variable t?  
    if tstop is None:
        info[3] = 0
    else:
        # check to make sure that tstop is greater than or equal to the final
        # time point requested
        if tstop < t[-1]:
            raise ValueError('tstop must be greater than or equal to the '
                             'final time point requested.')
        info[3] = 1
        rwork[0] = tstop
    # There appears to be a bug in ddaskr. If we're in intermediate output mode
    # and we have events, the daskr may skip a requested time and instead
    # return an event firing time. A reasonable fix appears to be using tstop
    # to prevent this.
    EVENT_INT_BUG_WORKAROUND = False
    if intermediate_output and nrt:
        EVENT_INT_BUG_WORKAROUND = True

    # Do you want the code to evaluate the partial derivatives automatically by
    # numerical differences?
    if jac is None:
        info[4] = 0
        jac = dummy_func
    else:
        info[4] = 1

    # info[5], info[6], info[7], info[8] are not used
    if hmax is not None:
        info[7] = 1
        rwork[2] = hmax

    # Do you want the code to solve the problem without invoking any special
    # inequality constraints?
    # LID is a variable used for setting up which variables are algebraic and
    # which are differential. The value of LID depends on info(9).
    # For now we only allow options 0 or 2 so that we don't have to worry
    # about telling the integrator how each variable is constrained.
    # Note:  This wrapper doesn't don't support info[9] = 0 or 2
    LID = 0
    if ineq_constr == False:
        info[9] = 0
        LID = 40
    #elif ineq_constr == 1:
    #    info[9] = 1
    #    LID = 40 + neq
    elif ineq_constr == True:
        info[9] = 2
        LID = 40
    #elif ineq_constr == 3:
    #    info[9] = 3
    #    LID = 40 + neq
    else:
        raise ValueError, 'ineq_constr must be True or False.'

    # Is it necessary to calculate the initial condition?
    # (i.e., are the initial conditions consistent?)
    if calculate_ic == False:
        info[10] = 0
    # if they are inconsistent, given Y_d, calculate Y_a and Y'_d must also
    # specify which components are algebraic and which are differential
    elif calculate_ic == True:
        info[10] = 1
        # check the var_types parameter
        if var_types is None:
            raise ValueError('if calculate_ic = 1, the var_types array '
                             'must be passed.')
        elif len(var_types) != len(y0):
            raise ValueError('var_types must be of length len(y0)')
        else:
            # +1 is for differential variables
            # -1 is for algebraic variables
            for ii in range(len(var_types)):
                if var_types[ii] == +1:
                    iwork[LID+ii] = +1
                elif var_types[ii] == -1:
                    iwork[LID+ii] = -1
                else:
                    raise ValueError('the var_types array may only contain '
                                     'entries of +1 or -1.')
                    

    # info[11] and info[12] are not used
    
    # Do you want to proceed to the integration after the initial condition
    # calculation is done?
    # (for now we will assume that we will always stop to capture the
    # calculated initial condition)
    # Note that once the initial condition is calculated, info[13] must
    # be reset to 0 to allow the integration to proceed
    info[13] = 1

    # info[14] is not used

    # Do you wish to control errors locally on all the variables?
    # (for now we will assume the answer is yes)
    info[15] = 0

    # Are the intial condition heuristic controls to be given their default
    # values?
    # (for now we will assume yes)
    info[16] = 0

    #info[18] is not used.

    # if no root function was passed in, we should assign the dummy function
    if rt is None:
        rt = dummy_func
        
    # assign the dummy function to psol since we don't use it
    psol = dummy_func
    
    
    # Set up for first call
    idid = 1 # First call
    t0, tindex = t[0], 1
    tout, yout = t0, y0
    t_root, y_root, i_root = [], [], []

    # variables for storing the output.
    tout_l = []
    yout_l = []
    ypout_l= []

    # add the initial point to the output array if the initial condition is
    # consistent
    if calculate_ic == False:
        tout_l.append(t0)
        yout_l.append(y0)
        ypout_l.append(yp0)

    # set the current time to the initial time
    tcurrent = t0

    # step_count is used to keep track of how many times idid = -1 occurs in a
    # row
    step_count = 0

    if redir_output == True:
        redir.start()

    if hasattr(res, '_cpointer'):
        res = res._cpointer
    if hasattr(jac, '_cpointer'):
        jac = jac._cpointer
    if hasattr(rt, '_cpointer'):
        rt = rt._cpointer

    try:
        while tcurrent < t[-1]:
            # set the desired output time
            twanted = t[tindex]
            if EVENT_INT_BUG_WORKAROUND:
                info[3] = 1
                if tstop and tstop < twanted:
                    rwork[0] = tstop
                else:
                    rwork[0] = twanted

            # daskr fails if tcurrent == twanted, but it's perfectly sensible
            # to just duplicate that point...
            if not tcurrent == twanted:
                try:
                    # continue the integration
                    treached, y, yp, idid, jroot = \
                            _daskr.ddaskr(res, tcurrent, y, yp, twanted,
                                          info, rtol, atol,
                                          rwork, iwork,
                                          rpar, ipar,
                                          jac, psol, rt, nrt)
                except (ValueError,OverflowError,TypeError,AttributeError,
                    AssertionError,FloatingPointError,ZeroDivisionError), e:
                    messages = redir.stop()
                    # report the error message
                    if messages is not None:
                        logger.warn(messages)
                    errMessage = "DASKR Error: Unrecoverable error: " + str(e)
                    logger.warn(errMessage)
                    outputs = (scipy.array(yout_l), scipy.array(tout_l),
                         scipy.array(ypout_l),
                         t_root, y_root, i_root)
                    raise daeintException(errMessage, outputs)

            # check for a negative value of idid so we know if there was a
            # problem
            if idid <= 0:
                # set appropriate options
                info[0]=1
                # idid=-1 is handled below
                if idid < -1:
                    messages = redir.stop()
                    # The task was interrupted
                    if messages is not None:
                        logger.warn(messages)
                    logger.warn(_msgs[idid])
                # if idid=-1, 500 steps have been taken sine the last time
                # idid=-1.
                # Whether we should continue depends on the value of max_steps.
                elif idid == -1:
                    # if we haven't hit the value of max_steps, then we should
                    # continue the integration.  Otherwise we should raise an
                    # exception.
                    if max_steps > step_count:
                        step_count += 500
                        continue
                    elif max_steps <= step_count:
                        messages = redir.stop()
                        # report the error message for idid = -1
                        if messages is not None:
                            logger.warn(messages)
                        logger.warn(_msgs[idid])

                   
                # send what output was obtained
                # If the integrator tried and failed to calculate the initial 
                #  condition, then yout_l will be of length 0, and have the
                #  wrong shape.
                # We replace it and ypout_l with arrays containing 0 rows, but
                #  the proper number of columns, for consistency with
                #  what we would get if it had succeeded.
                if len(yout_l) == 0:
                    yout_l = scipy.zeros((0, len(y0)), scipy.float_)
                    ypout_l = scipy.zeros((0, len(y0)), scipy.float_)
                outputs = (scipy.array(yout_l), scipy.array(tout_l),
                           scipy.array(ypout_l),
                           t_root, y_root, i_root)
                raise daeintException(_msgs[idid], outputs)
            
            # if there were no errors, continue to post output or take
            # appropriate action if we triggered an event or tstop
            else:

                # any time a succesful idid occurs we should reset step_count
                step_count = 0
                
                # updating of output should be done whether or not we are in
                # intermediate output mode

                # update the output
                tout_l.append(treached)
                yout_l.append(y)
                ypout_l.append(yp)

                # update the current time
                tcurrent = treached

                # if we hit one of the specified times, update time index
                if treached == t[tindex]:
                    # update the time index
                    tindex += 1

                # if the time we reached is tstop, then stop the integration
                if idid == 2:
                    break

                # if an initial condition was calculated, update info[13]
                # this allows the integration to proceed.  Also set
                # calculate_ic to False to avoid a redundant initial condition
                # calculation
                if idid == 4:
                    info[13] = 0
                    calculate_ic = False

                # if the solution was successful because a root of R(T,Y,Y') was
                # found at treached, then restart the integration.
                elif idid == 5:
                    t_root = tcurrent
                    y_root = y
                    i_root = jroot
                    # include the following break if we want all events to be
                    # terminating
                    break

                # Terminate if maximum number of timepoints has been reached.
                if len(tout_l) > max_timepoints:
                    messages = redir.stop()
                    # report the error message
                    if messages is not None:
                        logger.warn(messages)
                    errMessage = "DASKR Error: max_timepoints ("+               \
                        str(int(max_timepoints))+") exceeded."
                    logger.warn(errMessage)
                    outputs = (scipy.array(yout_l), scipy.array(tout_l),
                               scipy.array(ypout_l),
                               t_root, y_root, i_root)
                    raise daeintException(errMessage, outputs)

                    
        # format the output
        tout_l = scipy.array(tout_l)
        yout_l = scipy.array(yout_l)
        ypout_l = scipy.array(ypout_l)

        # Process outputs
        outputs = (yout_l, tout_l, ypout_l, t_root, y_root, i_root)

        return outputs

    finally:
        messages = redir.stop()
