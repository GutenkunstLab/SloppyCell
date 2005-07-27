_TEMP_DIR = '.SloppyCell'

import SloppyCell.Collections
KeyedList = SloppyCell.Collections.KeyedList
ExperimentCollection = SloppyCell.Collections.ExperimentCollection
Experiment = SloppyCell.Collections.Experiment
CalculationCollection = SloppyCell.Collections.CalculationCollection

import Ensembles
import Plotting
import Residuals
import Optimization
import Utility

# This bit of voodoo is to expose the Model class while also allowing the
#  module it's a part of to be reloaded. This preserves the convenience
#  of reload() for debugging while also making it easier to expose the Model
#  class.
#
# Note that code that wants to import SloppyCell.Model and get the module still
#  breaks.
if not locals().has_key('Model_mod'):
    import Model as Model_mod
Model = Model_mod.Model
