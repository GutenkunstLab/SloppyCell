_VERSION = 'CVS'
_TEMP_DIR = '.SloppyCell'

import logging
# I think I need to send one message to the root logger to get it started so the
#  lower level loggers can work.
logging.getLogger().debug('Logger started.')

# This is a slightly complicated, but nice, way to parse for a debug option.
from optparse import OptionParser
parser = OptionParser(version="SloppyCell version: %s" % _VERSION)
parser.add_option("--debug", dest="debug", metavar="FILE",
                  help="write debugging information to FILE. "
                  "If FILE is 'console' info will be sent to stderr")
(options, args) = parser.parse_args()
if options.debug:
    import Utility
    Utility.enable_debugging_msgs(options.debug)

import os
if not os.path.isdir(_TEMP_DIR): 
    os.mkdir(_TEMP_DIR)

#import Collections
#ExperimentCollection = Collections.ExperimentCollection
#Experiment = Collections.Experiment
#CalculationCollection = Collections.CalculationCollection
#
#import Ensembles
#import Plotting
#import Residuals
#import Observers
#import Optimization
#import Utility
#
## This bit of voodoo is to expose the Model and KeyedList classes while also 
##  allowing the modules they're in to be reloaded. This preserves the 
##  convenience of reload() for debugging while also making it easier to expose 
##  the classes.
##
## Note that code that wants to import SloppyCell.Model and get the module still
##  breaks.
#classimps = ['Model', 'KeyedList']
#for modname in classimps:
#    if not locals().has_key(modname+'_mod'):
#        exec 'import %s as %s_mod' % (modname, modname)
#    exec '%s = %s_mod.%s' % (modname, modname, modname)
