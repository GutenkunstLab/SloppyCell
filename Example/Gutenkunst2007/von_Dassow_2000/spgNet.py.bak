from SloppyCell.ReactionNetworks import *

# This file was distributed with Ingenue
f = file('spg1_4cell.net')
lines = f.readlines()


for ii, line in enumerate(lines):
    if line.startswith('&width'):
        width = int(line.split()[1])
    elif line.startswith('&height'):
        height = int(line.split()[1])
    elif line.startswith('&Network'):
        net_id = line.split()[1]
    elif line.startswith('&Genes'):
        begin_genes = ii
    elif line.startswith('&Interactions'):
        begin_ints = ii
    elif line.startswith('&ParameterValues'):
        begin_params = ii
    elif line.startswith('&InitLevels'):
        begin_inits = ii


net = Network(net_id)
if height > 1:
    for ii in range(height):
        for jj in range(width):
            net.add_compartment('cell_%i_%i' % (ii, jj))
else:
    for jj in range(width):
        net.add_compartment('cell_%i' % (jj))


ii = begin_genes + 1
while True:
    line = lines[ii].strip()
    if line.startswith('&endGenes'):
        break

    if line.startswith('&'):
        # Process this species
        species_id = line[1:]
        ii += 1
        # Skip to the end of this species entry
        while not lines[ii].strip().startswith('&end'):
            line = lines[ii]
            first, second = line.split()
            if first == '&Location':
                on_membrane = (second == 'membrane')
            ii += 1

        for comp_ii, comp_id in enumerate(net.compartments.keys()):
            if not on_membrane:
                net.add_species('%s_%s' % (species_id, comp_id), comp_id,
                                name = r'%s_{%i}' % (species_id, comp_ii))
            else:
                for jj in range(6):
                    id = '%s_%s_side_%i' % (species_id, comp_id, jj)
                    name = r'%s_{%i, %i}' % (species_id, comp_ii, jj)
                    net.add_species(id, comp_id, name = name)
    ii += 1

for comp_id in net.compartments.keys():
    net.set_var_constant('B_%s' % comp_id, True)

ii = begin_params + 1
while True:
    line = lines[ii].strip()
    if line.startswith('&endParameterValues'):
        break

    if line.startswith('&'):
        # Process this parameter
        temp = line.split()
        param_id, param_val = temp[0][1:], float(temp[1])
        net.add_parameter(param_id, param_val)
    ii += 1

# Create all the appropriate parameter names
for param_id in net.parameters.keys():
    if param_id == 'K_PTC_HH':
        name = r'k_{PTCHH}'
    elif param_id.startswith('K_'):
        term = param_id.split('_')[1]
        name = r'\kappa_{%s}' % term
    elif param_id.startswith('nu_'):
        term = param_id.split('_')[1]
        name = r'\nu_{%s}' % term
    elif param_id.startswith('H_'):
        term = param_id.split('_')[1]
        name = r'H_{%s}' % term
    elif param_id.startswith('alpha_'):
        term = param_id.split('_')[1]
        name = r'\alpha_{%s}' % term
    else:
        name = param_id
    net.parameters.get(param_id).name = name 

net.parameters.get('Endo_WG').name = r'r_{EndoWG}'
net.parameters.get('Exo_WG').name = r'r_{ExoWG}'
net.parameters.get('Mxfer_WG').name = r'r_{MxferWG}'
net.parameters.get('LMxfer_WG').name = r'r_{LMxferWG}'
net.parameters.get('LMxfer_PTC').name = r'r_{LMxferPTC}'
net.parameters.get('LMxfer_HH').name = r'r_{LMxferHH}'
net.parameters.get('maxHH').name = r'\left[HH\right]_0'
net.parameters.get('maxPTC').name = r'\left[PTC\right]_0'
net.parameters.get('C_CID').name = r'C_{CID}'
net.add_parameter('T_0', 1.0, is_optimizable=False)

ii = begin_inits + 1
while True:
    line = lines[ii].strip()
    if line.startswith('&endInitLevels'):
        break

    elif line.startswith('&BackgroundLevel'):
        spec_id = line.split()[1]
        value = float(line.split()[2])
        for var_id in net.species.keys():
            if var_id.startswith(spec_id):
                net.set_var_ic(var_id, value)

    elif line.startswith('&ColumnIC'):
        spec_id = lines[ii + 1].split()[1]
        value = float(lines[ii + 2].split()[1])
        column = int(lines[ii + 3].split()[1])
        cell_id = net.compartments.keys()[column]
        for var_id in net.species.keys():
            if var_id.startswith(spec_id) and var_id.count(cell_id):
                net.set_var_ic(var_id, value)
        ii += 3
    ii += 1

net.add_func_def('phi', ['X', 'k_X', 'v_X'], '(X**2)**(v_X/2) / ((k_X**2)**(v_X/2) + (X**2)**(v_X/2))',
                 name = r'\phi')
net.add_func_def('psi', ['X', 'k_X', 'v_X'], '1 - phi(X, k_X, v_X)',
                 name = r'\psi')

def presented_by_neighbors(net, comp_id, spec_id):
    ii = net.compartments.keys().index(comp_id)
    next = net.compartments.keys()[(ii + 1) % len(net.compartments)]
    prev = net.compartments.keys()[(ii - 1) % len(net.compartments)]
    terms = ['%s_%s_side_3' % (spec_id, comp_id),
             '%s_%s_side_4' % (spec_id, next),
             '%s_%s_side_5' % (spec_id, next),
             '%s_%s_side_0' % (spec_id, comp_id),
             '%s_%s_side_1' % (spec_id, prev),
             '%s_%s_side_2' % (spec_id, prev),
             ]
    return '+'.join(terms)

def total_in_cell(net, comp_id, spec_id):
    terms = ['%s_%s_side_%i' % (spec_id, comp_id, ii) for ii in range(6)]
    return '+'.join(terms)

def opposite_side(net, comp_id, side_jj, spec_id):
    ii = net.compartments.keys().index(comp_id)
    next = net.compartments.keys()[(ii + 1) % len(net.compartments)]
    prev = net.compartments.keys()[(ii - 1) % len(net.compartments)]
    if side_jj in [1, 2]:
        cell_id = next
    elif side_jj in [0, 3]:
        cell_id = comp_id
    elif side_jj in [4, 5]:
        cell_id = prev

    return '%s_%s_side_%i' % (spec_id, cell_id, (side_jj + 3)%6)


for comp_ii, comp_id in enumerate(net.compartments.keys()):
    # rhs for en_i
    #  First define EWG^{tot}_{n(i,j)}
    net.add_parameter('EWG_tot_pres_%s' % comp_id, is_optimizable=False, 
                      name=r'{EWG_{n(%i, j)}}^{tot}' % comp_ii)
    net.add_assignment_rule('EWG_tot_pres_%s' % comp_id, 
                            presented_by_neighbors(net, comp_id, 'EWG'))
    rule_str = 'T_0/H_en * (phi(EWG_tot_pres_%(comp)s * psi(CN_%(comp)s, K_CNen, nu_CNen), K_WGen, nu_WGen) - en_%(comp)s)' 
    net.add_rate_rule('en_%s' % comp_id, rule_str % {'comp': comp_id})

    # rhs for EN_i
    rule_str = 'T_0/H_EN * (en_%(comp)s - EN_%(comp)s)'
    net.add_rate_rule('EN_%s' % comp_id, rule_str % {'comp': comp_id})

    # rhs for wg_i
    num = '(beta_wg * phi(CID_%(comp)s * psi (CN_%(comp)s, K_CNwg, nu_CNwg), K_CIDwg, nu_CIDwg) + alpha_wg * phi(IWG_%(comp)s, K_WGwg, nu_WGwg))' % {'comp': comp_id}
    denom = '(1 + beta_wg * phi(CID_%(comp)s * psi(CN_%(comp)s, K_CNwg, nu_CNwg), K_CIDwg, nu_CIDwg) + alpha_wg * phi(IWG_%(comp)s, K_WGwg, nu_WGwg))' % {'comp': comp_id}
    rule_dict = {'num' : num, 'denom' : denom, 'comp' : comp_id}
    rule_str = 'T_0/H_wg * %(num)s/%(denom)s - T_0/H_wg * wg_%(comp)s'
    net.add_rate_rule('wg_%s' % comp_id, rule_str % rule_dict)

    # rhs for IWG_i
    net.add_parameter('EWG_tot_%s' % comp_id, is_optimizable=False, 
                      name=r'{EWG_{%i}}^{tot}' % comp_ii)
    net.add_assignment_rule('EWG_tot_%s' % comp_id, 
                            total_in_cell(net, comp_id, 'EWG'))
    rule_str = 'T_0/H_IWG * (wg_%(comp)s - IWG_%(comp)s) + T_0*(Endo_WG * EWG_tot_%(comp)s - Exo_WG * IWG_%(comp)s)'
    net.add_rate_rule('IWG_%s' % comp_id, rule_str % {'comp': comp_id})

    #rhs for EWG_i_j
    for side_jj in range(6):
        terms = ['T_0 * Exo_WG * IWG_%(comp)s/6',
                 '-(T_0 * Endo_WG * EWG_%(comp)s_side_%(side)i)',
                 'T_0 * Mxfer_WG * (%(sub)s - EWG_%(comp)s_side_%(side)i)',
                 'T_0 * LMxfer_WG * (EWG_%(comp)s_side_%(prev)i + EWG_%(comp)s_side_%(next)i - 2*EWG_%(comp)s_side_%(side)i)',
                 '-T_0/H_EWG * EWG_%(comp)s_side_%(side)i']
        rule_str = ' + '.join(terms)
        rule_dict =  {'sub': opposite_side(net, comp_id, side_jj, 'EWG'),
                      'comp': comp_id,
                      'side': side_jj,
                      'prev': (side_jj - 1) % 6,
                      'next': (side_jj + 1) % 6}
        net.add_rate_rule('EWG_%s_side_%i' % (comp_id, side_jj),
                          rule_str % rule_dict)

    # rhs for ptc_i
    rule_str = 'T_0/H_ptc * (phi(CID_%(comp)s * psi(CN_%(comp)s, K_CNptc, nu_CNptc), K_CIDptc, nu_CIDptc) - ptc_%(comp)s)'
    rule_dict = {'comp': comp_id}
    net.add_rate_rule('ptc_%s' % comp_id,
                      rule_str % rule_dict)

    # rhs for PTC_i_j
    for side_jj in range(6):
        terms = ['T_0/H_PTC * (ptc_%(comp)s/6 - PTC_%(comp)s_side_%(side)i)',
                 '-(T_0 * K_PTC_HH * maxHH * %(sub1)s * PTC_%(comp)s_side_%(side)i)',
                 'T_0 * LMxfer_PTC * (PTC_%(comp)s_side_%(prev)i + PTC_%(comp)s_side_%(next)i - 2*PTC_%(comp)s_side_%(side)i)']
        rule_str = ' + '.join(terms)
        rule_dict =  {'sub1': opposite_side(net, comp_id, side_jj, 'HH'),
                      'comp': comp_id,
                      'side': side_jj,
                      'prev': (side_jj - 1) % 6,
                      'next': (side_jj + 1) % 6}
        net.add_rate_rule('PTC_%s_side_%i' % (comp_id, side_jj),
                          rule_str % rule_dict)

    # rhs for cid_i
    rule_str = 'T_0/H_cid * (phi(B_%(comp)s * psi(EN_%(comp)s, K_ENcid, nu_ENcid), K_Bcid, nu_Bcid) - cid_%(comp)s)'
    rule_dict = {'comp': comp_id}
    net.add_rate_rule('cid_%s' % comp_id,
                      rule_str % rule_dict)

    # rhs for CID_i
    net.add_parameter('PTC_tot_%s' % comp_id, is_optimizable=False, 
                      name=r'{PTC_{%i}}^{tot}' % comp_ii)
    net.add_assignment_rule('PTC_tot_%s' % comp_id, 
                            total_in_cell(net, comp_id, 'PTC'))
    rule_str = 'T_0/H_CID * (cid_%(comp)s - CID_%(comp)s) - T_0 * C_CID * CID_%(comp)s * phi(PTC_tot_%(comp)s, K_PTCCID, nu_PTCCID)'
    rule_dict = {'comp': comp_id}
    net.add_rate_rule('CID_%s' % comp_id,
                      rule_str % rule_dict)

    # rhs for CN_i
    rule_str = 'T_0 * C_CID * CID_%(comp)s * phi(PTC_tot_%(comp)s, K_PTCCID, nu_PTCCID) - T_0 * CN_%(comp)s/H_CN'
    rule_dict = {'comp': comp_id}
    net.add_rate_rule('CN_%s' % comp_id,
                      rule_str % rule_dict)

    # rhs for hh_i
    rule_str = 'T_0/H_hh * (phi(EN_%(comp)s * psi(CN_%(comp)s, K_CNhh, nu_CNhh), K_ENhh, nu_ENhh) - hh_%(comp)s)'
    rule_dict = {'comp': comp_id}
    net.add_rate_rule('hh_%s' % comp_id,
                      rule_str % rule_dict)

    # rhs for HH_i_j
    for side_jj in range(6):
        terms = ['T_0/H_HH * (hh_%(comp)s/6 - HH_%(comp)s_side_%(side)s)'
                 '-T_0 * K_PTC_HH * maxPTC * %(sub1)s * HH_%(comp)s_side_%(side)s',
                 'T_0 * LMxfer_HH * (HH_%(comp)s_side_%(prev)i + HH_%(comp)s_side_%(next)i - 2*HH_%(comp)s_side_%(side)i)']
        rule_str = ' + '.join(terms)
        rule_dict =  {'sub1': opposite_side(net, comp_id, side_jj, 'PTC'),
                      'comp': comp_id,
                      'side': side_jj,
                      'prev': (side_jj - 1) % 6,
                      'next': (side_jj + 1) % 6}
        net.add_rate_rule('HH_%s_side_%i' % (comp_id, side_jj),
                          rule_str % rule_dict)

    # rhs for PH_i_j
    for side_jj in range(6):
        rule_str= 'T_0 * K_PTC_HH * maxHH * %(sub1)s * PTC_%(comp)s_side_%(side)i - T_0*PH_%(comp)s_side_%(side)i / H_PH'
        rule_dict =  {'sub1': opposite_side(net, comp_id, side_jj, 'HH'),
                      'comp': comp_id,
                      'side': side_jj}
        net.add_rate_rule('PH_%s_side_%i' % (comp_id, side_jj),
                          rule_str % rule_dict)
