import SloppyCell.Ensembles as Ensembles
import SloppyCell.Residuals as Residuals
import SloppyCell.Observers as Observers
import SloppyCell.Optimization as Optimization
import SloppyCell.Utility as Utility

import SloppyCell.KeyedList_mod as KeyedList_mod
KeyedList = KeyedList_mod.KeyedList
import SloppyCell.Model_mod as Model_mod
Model = Model_mod.Model
import SloppyCell.Collections as Collections
Experiment = Collections.Experiment

try:
    import Dynamics
except ImportError:
    pass

import IO

try:
    import Plotting
except ImportError:
    pass

import Reactions
import PerfectData

import Network_mod
Network = Network_mod.Network

# This is the same voodoo that's used in SloppyCell/__init__.py
#classimps = ['Network']
#for modname in classimps:
#    if not locals().has_key(modname+'_mod'):
#        exec 'import %s as %s_mod' % (modname, modname)
#    exec '%s = %s_mod.%s' % (modname, modname, modname)
