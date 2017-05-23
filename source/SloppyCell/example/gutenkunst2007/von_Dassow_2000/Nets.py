import scipy

from SloppyCell.ReactionNetworks import *

import spgNet
net = spgNet.net
net.compile()

# Load the parameter set from spg.params, which was included with Ingenue
f = file('spg.params')
lines = f.read().split('\r')
f.close()
param_sets = []
param_scores = []
lineiter = lines.__iter__()
for line in lineiter:
    if line.strip().startswith('&par'):
        params = KeyedList()
        line = lineiter.next().strip()
        while not line:
            line = lineiter.next().strip()


        # Skip the blank and Score info
        score_line = lineiter.next()
        words = score_line.strip().split()
        param_scores.append(float(words[1]))

        line = lineiter.next()
        line = lineiter.next()
        while not line.startswith('&endFPARS'):
            words = line.strip().split()
            params.set(words[0][1:], float(words[1]))
            line = lineiter.next()
        param_sets.append(params)
param_scores = scipy.asarray(param_scores)

# For our network we use the set of parameters that gave the lowest score
set_ii = scipy.argmin(param_scores)
net1 = net.copy('spg_%i' % set_ii)
net1.update_optimizable_vars(param_sets[set_ii])

networks = [net1]
int_times = [(0, 1000)] * len(networks)
