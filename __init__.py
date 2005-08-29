_TEMP_DIR = '.SloppyCell'

import Collections
ExperimentCollection = Collections.ExperimentCollection
Experiment = Collections.Experiment
CalculationCollection = Collections.CalculationCollection

import Ensembles
import Plotting
import Residuals
import Observers
import Optimization
import Utility

# This bit of voodoo is to expose the Model and KeyedList classes while also 
#  allowing the modules they're in to be reloaded. This preserves the 
#  convenience of reload() for debugging while also making it easier to expose 
#  the classes.
#
# Note that code that wants to import SloppyCell.Model and get the module still
#  breaks.
classimps = ['Model', 'KeyedList']
for modname in classimps:
    if not locals().has_key(modname+'_mod'):
        exec 'import %s as %s_mod' % (modname, modname)
    exec '%s = %s_mod.%s' % (modname, modname, modname)
