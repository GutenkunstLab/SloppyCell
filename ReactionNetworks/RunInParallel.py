from SloppyCell.RunInParallel import *

# Import the necessary items into the namespace for the worker.
statement_to_all_workers('from SloppyCell.ReactionNetworks.Network_mod import Network')
statement_to_all_workers('import SloppyCell.ReactionNetworks.Dynamics as Dynamics')
statement_to_all_workers('import SloppyCell.ReactionNetworks.PerfectData as PerfectData')
statement_to_all_workers('import SloppyCell.Ensembles as Ensembles')
