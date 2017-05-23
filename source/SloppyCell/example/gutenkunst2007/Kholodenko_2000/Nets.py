from SloppyCell.ReactionNetworks import *

net1 = IO.from_SBML_file('BIOMD0000000010.xml')
net1.set_var_optimizable('n', False)

net2 = net1.copy('Cooperative')
p2 = KeyedList([('n', 2),
                ('Ki',  18),
                ('K1', 50),
                ('KK2', 40),
                ('KK3', 100),
                ('KK4', 100),
                ('KK5', 100),
                ('KK6', 100),
                ('KK7', 100),
                ('KK8', 100),
                ('KK9', 100),
                ('KK10', 100),
                ('V9', 1.25),
                ('V10', 1.25),
                ])
for id, value in p2.items():
    net2.set_var_ic(id, value)

networks = [net1, net2]
int_times = [(0, 150*60), (0, 205*60)]
