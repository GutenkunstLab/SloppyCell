from SloppyCell.ReactionNetworks import *

# Load our network from an SBML file
base_net = IO.from_SBML_file('BIOMD0000000005.xml')

###############################################################################

# That SBML file is a little sloppy. Let's clean up a couple things.

# First, we force the 'cell', and 'EmptySet' variables to have constant values.
# This doesn't change the dynamics at all, but does speed up the integration
#  a little bit.
base_net.set_var_constant('cell', True)
base_net.set_var_constant('EmptySet', True)

# We'll define the 'total cyclin' species that is used in plots.
base_net.add_species('YT', 'cell', name = 'total cyclin')
base_net.addAssignmentRule('YT', 'Y+YP+pM+M')

# These two variables are set to 0 in the original model. Let's make sure they
#  stay that way when we optimize.
base_net.set_var_optimizable('k5notP', False)
base_net.set_var_optimizable('k2', False)

###############################################################################

# Now we'll copy our network so we can mess with it and not affect the original.
# The 'perturbed' argument gives this network a new id, so it can be used in
#  a Model without running over the base_network.
perturbed_net = base_net.copy('perturbed')

# We'll change one of our parameters.
perturbed_net.set_initial_var_value('k6', 2.0)
# Let's make sure k6 can't be optimized in this network.
perturbed_net.set_var_optimizable('k6', False)

# We want to start from a fixed point. 
#
# First we integrate for a while.
times = [0, 100]
traj = perturbed_net.integrate(times)
# Then we find the nearest fixed point
fp = perturbed_net.dyn_var_fixed_point()
# Finally, we set our dynamic variables initial conditions to that fixed point.
perturbed_net.set_dyn_var_ics(fp)

# We add a few perturbations to our network
perturbed_net.addEvent('perturb1', 'gt(time, 10)', {'M': '1.2*M'})
perturbed_net.addEvent('perturb2', 'gt(time, 30)', {'M': '2*M'})
perturbed_net.addEvent('perturb3', 'gt(time, 60)', {'M': '2.4*M'})

###############################################################################

# This network will simulate cell growth as done in the paper, by allowing k6
#  to decline exponentially, then double upon division.
growth_net = base_net.copy('growth')

# Add a division time parameter
growth_net.add_parameter('Td', 116)

# We'll make k6 decay exponentially
growth_net.add_rate_rule('k6', '-0.693 * k6/Td')
growth_net.set_initial_var_value('k6', 2.5)

# But if k6 gets less than 1.3, the cell divides and k6 doubles
growth_net.addEvent('gene_replication', 'lt(k6, 1.3)', {'k6': '2 * k6'})

networks = [base_net, perturbed_net, growth_net]
int_times = [(0,100), (0,100), (0, 500)]
