_VERSION = 'CVS'
import os
_TEMP_DIR = os.path.join(os.getcwd(), '.SloppyCell')

import logging
logging.basicConfig()
logger = logging.getLogger('__init__')

# Check for debugging option. I tried using optparse for this, but ran into
# issues with ipython and mpirun, both of which pollute sys.argv.
disable_c = False
import sys
args_to_remove = []
for arg in sys.argv:
    if arg.startswith('--debugSC'):
        words = arg.split('=')
        import Utility
        if len(words) == 2:
            Utility.enable_debugging_msgs(words[1])
        else:
            Utility.enable_debugging_msgs(None)
        args_to_remove.append(arg)
    elif arg.startswith('--disableC'):
        disable_c = True
        args_to_remove.append(arg)
currdir = os.getcwd()
# We need to remove these arguments from the list before they get to uniitest
#  or it will complain.
for arg in args_to_remove:
    sys.argv.remove(arg)
            
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    os.chdir(currdir)
    HAVE_MPI = True
    num_procs = comm.Get_size()
    my_rank = comm.Get_rank()
    my_host = MPI.Get_processor_name()
except ImportError:
    os.chdir(currdir)
    HAVE_MPI = False
    num_procs = 1
    my_rank = 0
    comm = None
    import socket
    my_host = socket.gethostname()
logger.debug('Node %i is on host %s.' % (my_rank, my_host))

if my_rank == 0 and not os.path.isdir(_TEMP_DIR):
    os.mkdir(_TEMP_DIR)
