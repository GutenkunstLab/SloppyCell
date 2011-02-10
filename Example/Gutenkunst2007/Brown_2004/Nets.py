from SloppyCell.ReactionNetworks import *

# Load the network from the XML file
net = IO.from_SBML_file('BIOMD0000000033.xml')
net.set_var_ic('EGF', 0)
net.set_var_ic('NGF', 0)
net.compile()

# Create our various conditions
# EGF stimulation in the wild-type
egf_stim = net.copy('EGF_stim')
egf_stim.set_var_ic('EGF', 10002000)

# NGF stimulation in the wild-type
ngf_stim = net.copy('NGF_stim')
ngf_stim.set_var_ic('NGF', 456000)

# EGF stimulation with PI3K inhibited
egf_ly = egf_stim.copy('EGF_LY')
egf_ly.set_var_ic('kPI3KRas', 0)
egf_ly.set_var_optimizable('kPI3KRas', False)
egf_ly.set_var_ic('kPI3K', 0)
egf_ly.set_var_optimizable('kPI3K', False)

# NGF stimulation with PI3K inhibited
ngf_ly = ngf_stim.copy('NGF_LY')
ngf_ly.set_var_ic('kPI3KRas', 0)
ngf_ly.set_var_optimizable('kPI3KRas', False)
ngf_ly.set_var_ic('kPI3K', 0)
ngf_ly.set_var_optimizable('kPI3K', False)

# NGF stimulation with a dominant-negative Rap1
ngf_DN_Rap1 = ngf_stim.copy('NGF_DN_Rap1')
ngf_DN_Rap1.set_var_ic('kRap1ToBRaf', 0)
ngf_DN_Rap1.set_var_optimizable('kRap1ToBRaf', False)

# NGF stimulation with a dominant-negative Ras
ngf_DN_Ras = ngf_stim.copy('NGF_DN_Ras')
ngf_DN_Ras.set_var_ic('kRasToRaf1', 0)
ngf_DN_Ras.set_var_optimizable('kRasToRaf1', False)
ngf_DN_Ras.set_var_ic('kPI3KRas', 0)
ngf_DN_Ras.set_var_optimizable('kPI3KRas', False)

# Collect our networks in a list
networks = [egf_stim, ngf_stim, egf_ly, ngf_ly, ngf_DN_Rap1, ngf_DN_Ras]
# And make a list of corresponding integration times
int_times = [[0, 120]] * len(networks)
