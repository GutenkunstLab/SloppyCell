from SloppyCell.ReactionNetworks import *

net = Network('frog_cell')
net.add_compartment('cell')

net.add_species('M', 'cell', 0.0)
net.add_rate_rule('M', '(v_dp*(1-D) + v_dpp*D)*(C_T - M) - (v_wp*(1-W) + v_wpp*W)*M')
net.add_species('D', 'cell', 0.0)
#net.add_rate_rule('D', 'v_d*(M*(1-D)/(K_md/mu + (1-D)) - rho_d*D/(K_mdr/mu + D))')
net.add_rate_rule('D', 'v_d*(M*(1-D) - rho_d*D/(K_mdr/mu + D))')

net.add_species('W', 'cell', 0.0)
net.add_rate_rule('W', 'v_w*(-M*W/(K_mw/mu + W) + rho_w*(1-W)/(K_mwr/mu + (1-W)))')

net.add_parameter('v_dp', 0.0152)
net.add_parameter('v_dpp', 0.187)
net.add_parameter('rho_d', 0.0079)
net.add_parameter('v_wp', 6.4e-7)
net.add_parameter('v_wpp', 1.05)
net.add_parameter('rho_w', 0.0366)
net.add_parameter('K_mdr', 0.111)
net.add_parameter('K_mw', 0.0196)
net.add_parameter('K_mwr', 0.0736)
net.add_parameter('v_d', 7.24)
net.add_parameter('v_w', 2.21)

net.add_parameter('mu', 1.0, is_optimizable=False)
net.add_parameter('C_T', is_optimizable=False)

A = net.copy('net_A')
A.add_species('L', 'cell', 0.0)
A.add_rate_rule('L', '(v_wp*(1-W) + v_wpp*W)*(1-L)')
A.set_var_ic('C_T', 0)
A.set_var_ic('W', 1)

B = net.copy('net_B')
B.add_species('L', 'cell', 0.0)
B.add_rate_rule('L', '(v_wp*(1-W) + v_wpp*W)*(1-L)')
B.set_var_ic('C_T', 1)
B.set_var_ic('M', 1)
B.set_var_ic('D', 1)
# I was getting some integration problems in the hessian calculation. They
#  seem to be caused by W immediately shooting to this value, so I just set the
#  IC there.
B.set_var_ic('W', 1.5e-3)

C = net.copy('net_C')
C.add_species('L_p', 'cell', 1.0)
C.add_rate_rule('L_p', '-(v_dp*(1-D) + v_dpp*D)*L_p')
C.set_var_ic('C_T', 0)
C.set_var_ic('W', 1)

D = net.copy('net_D')
D.add_species('L_p', 'cell', 1.0)
D.add_rate_rule('L_p', '-(v_dp*(1-D) + v_dpp*D)*L_p')
D.set_var_ic('C_T', 1)
D.set_var_ic('M', 1)
D.set_var_ic('D', 1)
D.set_var_ic('W', 1.5e-3)

E = net.copy('net_E')
E.set_var_ic('C_T', 0.25)
E.set_var_ic('M', 0.25)
E.set_var_ic('mu', 0.83)

F = net.copy('net_F')
F.set_var_ic('C_T', 0)
F.set_var_ic('D', 1.0)
F.set_var_ic('mu', 0.83)

G = net.copy('net_G')
G.set_var_ic('C_T', 0.254)
G.set_var_ic('M', 0.254)
G.set_var_ic('D', 1.0)
G.set_var_ic('W', 1.0)
G.set_var_ic('mu', 0.67)

H = net.copy('net_H')
H.set_var_ic('C_T', 0)
H.set_var_ic('D', 1.0)
H.set_var_ic('mu', 0.67)

networks = [A, B, C, D, E, F, G, H]
int_times = [(0, 20), (0, 20), (0, 20), (0, 20), (0, 20), (0, 40), (0, 20), 
             (0, 20)]
