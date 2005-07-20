from SloppyCell import *

import IO
import Plotting
import Reactions

# This is the same voodoo that's used in SloppyCell/__init__.py
if not locals().has_key('Network_mod'):
    import Network as Network_mod
Network = Network_mod.Network
