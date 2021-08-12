import scipy

from SloppyCell.ReactionNetworks import *

import LeeNet
reload(LeeNet)
from LeeNet import net

# Network that is at unstimulated fixed point
net.set_var_ic('W', 0)
traj1 = Dynamics.integrate(net, [0, 1e5], rtol=1e-12)
unstimulated = net.copy('unstimulated')
for var in unstimulated.dynamicVars.keys():
    if not isinstance(net.get_var_ic(var), str):
        unstimulated.set_var_ic(var, traj1.get_var_val_index(var, -1))

fig2a = unstimulated.copy('fig2a')
fig2a.set_var_ic('v12', 0)
fig2a.set_var_optimizable('v12', False)
fig2a.set_var_ic('v14', 0)
fig2a.set_var_optimizable('v14', False)

fig2b = unstimulated.copy('fig2b')
fig2b.set_var_ic('v12', 0)
fig2b.set_var_optimizable('v12', False)
fig2b.set_var_ic('v14', 0)
fig2b.set_var_optimizable('v14', False)
fig2b.set_var_ic('X12', 0.2)

fig2c = unstimulated.copy('fig2c')
fig2c.set_var_ic('v12', 0)
fig2c.set_var_optimizable('v12', False)
fig2c.set_var_ic('v14', 0)
fig2c.set_var_optimizable('v14', False)
fig2c.set_var_ic('k2', 0)
fig2c.set_var_optimizable('k2', False)
fig2c.set_var_ic('X2', 1000)
# We also need to make sure Dsh0 is updated.
fig2c.set_var_ic('Dsh0', 1100)

fig2d = unstimulated.copy('fig2d')
fig2d.set_var_ic('v12', 0)
fig2d.set_var_optimizable('v12', False)
fig2d.set_var_ic('v14', 0)
fig2d.set_var_optimizable('v14', False)
fig2d.set_var_ic('k4', 0)
fig2d.set_var_optimizable('k4', False)
fig2d.set_var_ic('k9', 0)
fig2d.set_var_optimizable('k9', False)

fig2e = unstimulated.copy('fig2e')
fig2e.set_var_ic('v12', 0)
fig2e.set_var_optimizable('v12', False)
fig2e.set_var_ic('v14', 0)
fig2e.set_var_optimizable('v14', False)
fig2e.set_var_ic('TCF0', 1000)
# X11, X13, and X14 are assumed to be in equilibrium. Adding TCF changes this
#  equilibrium and thus changes the initial concentration of X11 we need to
#  use for this integration.
# This number results from constraining the total amount of non-active BCatenin
#  to be equal to the steady state unstimulated value: (Denoted BcNA)
BcNA = unstimulated.get_var_ic('X11')
TCF0 = fig2e.get_var_ic('TCF0')
K16 = fig2e.get_var_ic('K16')
val = 0.5*(BcNA - K16 - TCF0 + scipy.sqrt(4*BcNA*K16 + (-BcNA+K16+TCF0)**2))
fig2e.set_var_ic('X11', val)

fig6a = unstimulated.copy('transient_a')
fig6a.add_parameter('lam', 1./20, is_optimizable=False)
fig6a.add_assignment_rule('W', 'exp(-lam*time)')

fig6b = fig6a.copy('transient_b')
fig6b.set_var_ic('v14', fig6a.get_var_ic('v14')*5)
fig6b.set_var_ic('k15', fig6a.get_var_ic('k15')*5)

fig6c = fig6a.copy('transient_c')
fig6c.set_var_ic('v14', fig6a.get_var_ic('v14')/5)
fig6c.set_var_ic('k15', fig6a.get_var_ic('k15')/5)


networks = [fig2a, fig2b, fig2c, fig2d, fig2e, fig6a, fig6b, fig6c]
int_times = [(0, 3*60), (0, 3*60), (0, 3*60), (0, 3*60), (0, 3*60), 
             (0, 16*60), (0, 16*60), (0, 16*60)]
