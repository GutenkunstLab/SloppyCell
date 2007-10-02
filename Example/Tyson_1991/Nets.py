from SloppyCell.ReactionNetworks import *

# Load our network from an SBML file
base_net = IO.from_SBML_file('BIOMD0000000005.xml', 'base')

###############################################################################

# That SBML file is a little sloppy. Let's clean up a couple things.

# First, we force the 'cell', and 'EmptySet' variables to have constant values.
# This doesn't change the dynamics at all, but does speed up the integration
#  a little bit.
base_net.set_var_constant('cell', True)
base_net.set_var_constant('EmptySet', True)

# We'll define the 'total cyclin' species that is used in plots.
base_net.add_species('YT', 'cell', name = 'total cyclin')
base_net.add_assignment_rule('YT', 'Y+YP+pM+M')

# These two variables are set to 0 in the original model. Let's make sure they
#  stay that way when we optimize.
base_net.set_var_optimizable('k5notP', False)
base_net.set_var_optimizable('k2', False)

###############################################################################

# Output the LaTeX'd equations. This can be very useful for debugging.
IO.eqns_TeX_file(base_net, 'eqns.tex')

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
Dynamics.integrate(perturbed_net, [0, 1000])
# Then we find the nearest fixed point
fp = Dynamics.dyn_var_fixed_point(perturbed_net)
# Finally, we set our dynamic variables initial conditions to that fixed point.
perturbed_net.set_dyn_var_ics(fp)

# We add a few perturbations to our network
perturbed_net.add_event('perturb1', 'gt(time, 10)', {'M': '1.2*M'})
perturbed_net.add_event('perturb2', 'gt(time, 30)', {'M': '2*M'})
perturbed_net.add_event('perturb3', 'gt(time, 60)', {'M': '2.4*M'})

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
growth_net.add_event('gene_replication', 'lt(k6, 1.3)', {'k6': '2 * k6'})

###############################################################################

my_net = base_net.copy('mine')

# We add a new species to the network called killer, with an initial 
#  concentration of 2.0
my_net.add_species('killer', 'cell', 2.0)

# Let's make killer dimerize irreversibly with pM
my_net.add_species('killpM', 'cell', 0.0)
my_net.add_parameter('kf_kpM', 1.0)
my_net.addReaction(Reactions.HeterodimerizationReaction,
                    'kill_pM_binding',
                    dimer='killpM',
                    A = 'pM', B = 'killer',
                    rate = 'kf_kpM')

# Free killer also degrades in an unusual time-dependent manner
my_net.addReaction('killer_degredation', {'killer': -1},
                    '0.1*killer * sin(time/25)')
