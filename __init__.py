_VERSION = 'CVS'
_TEMP_DIR = '.SloppyCell'

import logging
logging.basicConfig()

# Check for debugging option. I tried using optparse for this, but ran into
# issues with ipython and mpirun, both of which pollute sys.argv.
import sys
for arg in sys.argv:
    if arg.startswith('--debugSC'):
        import Utility
        words = arg.split('=')
        if len(words) == 2:
            Utility.enable_debugging_msgs(words[1])
        else:
            Utility.enable_debugging_msgs(None)
            
import os
if not os.path.isdir(_TEMP_DIR): 
    os.mkdir(_TEMP_DIR)
