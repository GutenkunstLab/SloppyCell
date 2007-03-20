import logging
logger = logging.getLogger('Parallel')

import sys, traceback

import SloppyCell
import Collections

from SloppyCell import num_procs, my_rank, my_host, HAVE_PYPAR
if HAVE_PYPAR:
    import pypar

import SloppyCell.Utility as Utility

class Statement:
    """
    Class for sending Python statements to workers.
    """
    def __init__(self, statement, locals={}):
        self.statement = statement
        self.locals = locals

while my_rank != 0:
    # Wait for a message
    message = pypar.receive(source=0)

    # If the message is a SystemExit exception, exit the code.
    if isinstance(message, SystemExit):
        sys.exit()

    # Exception handling:
    #    If we catch a SloppyCellException during a eval(), it's probably just 
    #      a numerical issue so we just pass it back to the master to deal with.
    #      Note that we don't catch SloppyCellExceptions for exec'd things.
    #      This is because exec'd things shouldn't return anything, thus the
    #      master won't be waiting for a reply.
    #    If we catch any other exception, it's probably a bug in the code. Print
    #      a nice traceback, save results, and exit the code.
    try:
        if isinstance(message, Statement):
            command, msg_locals = message.statement, message.locals
            locals().update(msg_locals)
            exec(command)
        else:
            command, msg_locals = message 
            locals().update(msg_locals)
            try:
                result = eval(command)
                pypar.send(result, 0)
            except Utility.SloppyCellException, X:
                pypar.send(X, 0)
    except:
        # Assemble and print a nice traceback
        tb = traceback.format_exception(sys.exc_type, sys.exc_value, 
                                        sys.exc_traceback)
        logger.critical(('node %i:'%my_rank).join(tb))
        save_to = '.SloppyCell/node_%i_crash.bp' % my_rank
        logger.critical("node %i: Command being run was: %s."
                        % (my_rank, command))
        Utility.save(msg_locals, save_to)
        logger.critical("node %i: Corresponding locals saved to %s."
                        % (my_rank, save_to))
        sys.exit()

def stop_workers():
    """
    Send all workers the command to exit the program.
    """
    for worker in range(1, num_procs):
        pypar.send(SystemExit(), worker)

if my_rank == 0:
    import atexit
    atexit.register(stop_workers)

def statement_to_all_workers(statement, locals={}):
    """
    Send a Python statement to all workers for execution.
    """
    for worker in range(1, num_procs):
        pypar.send(Statement(statement, locals), worker)
