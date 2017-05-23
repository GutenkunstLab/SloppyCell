from SloppyCell.ReactionNetworks import *

heme_net = Network('hemeglobin')
heme_net.add_compartment('test_tube')
heme_net.add_parameter('c', 0, name='free_ligand', typical_value=100)
heme_net.add_rate_rule('c', '1')
heme_net.add_parameter('K1', 0.0163)
heme_net.add_parameter('K2', 0.0556)
heme_net.add_parameter('K3', 0.00370)
heme_net.add_parameter('K4', 0.499)
heme_net.add_species('ligand_per_protein', 'test_tube', typical_value=2)
heme_net.add_assignment_rule('ligand_per_protein', '(c*K1 + 2*c**2*K1*K2 + 3*c**3*K1*K2*K3 + 4*c**4*K1*K2*K3*K4)/(1 + c*K1 + c**2*K1*K2 + c**3*K1*K2*K3 + c**4*K1*K2*K3*K4)')

alb_net = Network('albumin')
alb_net.add_compartment('test_tube')
alb_net.add_parameter('c', 0, name='free_ligand', typical_value=100)
alb_net.add_rate_rule('c', '1')
alb_net.add_parameter('K1', 17716)
alb_net.add_parameter('K2', 16489)
alb_net.add_parameter('K3', 1866)
alb_net.add_parameter('K4', 2940)
alb_net.add_parameter('K5', 117)
alb_net.add_parameter('K6', 957)
alb_net.add_species('ligand_per_protein', 'test_tube', typical_value=3)
alb_net.add_assignment_rule('ligand_per_protein', '(c*K1 + 2*c**2*K1*K2 + 3*c**3*K1*K2*K3 + 4*c**4*K1*K2*K3*K4 + 5*c**5*K1*K2*K3*K4*K5 + 6*c**6*K1*K2*K3*K4*K5*K6)/(1 + c*K1 + c**2*K1*K2 + c**3*K1*K2*K3 + c**4*K1*K2*K3*K4 + c**5*K1*K2*K3*K4*K5 + c**6*K1*K2*K3*K4*K5*K6)')

networks = [alb_net, heme_net]
int_times = [(1e-6, 10000e-6), (0.1, 1000)]
