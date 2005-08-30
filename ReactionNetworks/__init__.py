from SloppyCell import *

try:
    import Dynamics
except:
    pass
import IO
import Plotting
import Reactions
import PerfectData

# This is the same voodoo that's used in SloppyCell/__init__.py
classimps = ['Network']
for modname in classimps:
    if not locals().has_key(modname+'_mod'):
        exec 'import %s as %s_mod' % (modname, modname)
    exec '%s = %s_mod.%s' % (modname, modname, modname)
