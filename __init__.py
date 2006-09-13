_VERSION = 'CVS'
_TEMP_DIR = '.SloppyCell'

import logging
logging.basicConfig()
logger = logging.getLogger('__init__')

# Check for debugging option. I tried using optparse for this, but ran into
# issues with ipython and mpirun, both of which pollute sys.argv.
import sys
for arg in sys.argv:
    if arg.startswith('--debugSC'):
        words = arg.split('=')
        import Utility
        if len(words) == 2:
            Utility.enable_debugging_msgs(words[1])
        else:
            Utility.enable_debugging_msgs(None)
            
try:
    import pypar
    HAVE_PYPAR = True
    num_procs = pypar.size()
    my_rank = pypar.rank()
    my_host = pypar.get_processor_name()
    import atexit
    atexit.register(pypar.finalize)
except:
    HAVE_PYPAR = False
    num_procs = 1
    my_rank = 0
    import socket
    my_host = socket.gethostname()
logger.debug('Node %i is on host %s.' % (my_rank, my_host))

import os
if my_rank == 0 and not os.path.isdir(_TEMP_DIR): 
    os.mkdir(_TEMP_DIR)

import OldScipySupport
