from SloppyCell.ReactionNetworks import *

net = IO.from_SBML_file('BIOMD0000000035.xml', 'base')

base = net.copy('base')

# For figure 5
stable = net.copy('stable')
stable.set_var_ic('deltaR', 0.05)

networks = [net, stable]
int_times = [(0, 50), (0, 100)]

reduced = Network('reduced')
reduced.add_compartment('base')
for id in net.parameters.keys():
    reduced.add_parameter(id, net.get_var_ic(id))

reduced.add_species('R', 'base', 0)
reduced.add_species('C', 'base', 0)

reduced.add_parameter('Kd', is_constant=False)
reduced.add_assignment_rule('Kd', 'thetaA/gammaA')

reduced.add_parameter('rho', is_constant=False)
reduced.add_assignment_rule('rho', 'betaA/(deltaMA * (gammaC*R + deltaA))')

reduced.add_parameter('A_tilde', is_constant=False)
reduced.add_assignment_rule('A_tilde', '.5*(alphaAp * rho - Kd) + .5 * sqrt((alphaAp * rho - Kd)**2 + 4*alphaA*rho*Kd)')

reduced.add_rate_rule('R', 'betaR/deltaMR * (alphaR*thetaR + alphaRp*gammaR*A_tilde)/(thetaR + gammaR*A_tilde) - gammaC*A_tilde*R + deltaA*C - deltaR*R')
reduced.add_rate_rule('C', 'gammaC*A_tilde*R - deltaA * C')
