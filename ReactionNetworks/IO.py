import os

try:
    import SBMLInterface
    from SBMLInterface import toSBMLFile as to_SBML_file
    from SBMLInterface import fromSBMLFile as from_SBML_file
except ImportError:
    print 'SBML import and export not available.'

import SloppyCell
import SloppyCell.ExprManip as ExprManip
import Network_mod

def net_DOT_file(net, filename = None):
    lines = []
    lines.append('digraph "%s" {' % net.id)
    lines.append('\tsize="7.5,10!"')
    for id in net.species.keys():
        lines.append('\t"%s"[color=black]' % net.get_component_name(id))

    lines.append('')
    for rid, rxn in net.reactions.items():
        rxn_name = net.get_component_name(rid)
        lines.append('\t"%s"[shape=rectangle][color=red]' % rxn_name)
        for rid, stoich in rxn.stoichiometry.items():
            rname = net.get_component_name(rid)
            if stoich < 0:
                lines.append('\t\t"%s" -> "%s";' % (rxn_name, rname))
            elif stoich > 0:
                lines.append('\t\t"%s" -> "%s";' % (rname, rxn_name))
            else:
                lines.append('\t\t"%s" -> "%s"[arrowhead=dot];' % (rname, 
                                                                   rxn_name))
        
    lines.append('}')

    if filename is None:
        filename = '%s.dot' % net.id
    f = file(filename, 'w')
    f.write(os.linesep.join(lines))
    f.close()


def eqns_TeX_file(net, filename = None, sub_funcdefs = True):
    net.compile()

    lines = []
    lines.append(r'\documentclass{article}')
    lines.append(r'\usepackage{amsmath}')
    lines.append(r'\usepackage{fullpage}')
    lines.append(r'\usepackage{longtable}')
    lines.append(r'\begin{document}')
    lines.append(_net_eqns_to_TeX(net, sub_funcdefs))
    lines.append(r'\end{document}')

    if filename is None:
        filename = '%s.tex' % net.id
    f = file(filename, 'w')
    f.write(os.linesep.join(lines))
    f.close()

def _net_eqns_to_TeX(net, sub_funcdefs):
    """
    Return a string that contains the longtable-bound TeX'd equations for the network
    """
    # Build up our name_dict. We wrap variables in a mathrm
    name_dict = dict([(id, r'\mathrm{%s}' % net.get_component_name(id, True))
                      for id in net.variables.keys()] + 
                     [(id, net.get_component_name(id, True))
                       for id in net.functionDefinitions.keys()])
    # Species get wrapped in [ ]
    species_dict = dict([(id, r'\left[\mathrm{%s}\right]'
                          % net.get_component_name(id, True))
                         for id in net.species.keys()])
    name_dict.update(species_dict)
    eqns_str = ExprManip.Py2TeX.dict2TeX(net.diff_eq_rhs, name_dict, 
                                       r'\frac{d\,%s}{dt}', longtable=True)
    assigns_str = ExprManip.Py2TeX.dict2TeX(net.assignmentRules, name_dict, 
                                          longtable=True)

    return os.linesep.join([eqns_str, assigns_str])

def dynamic_function_from_file(obj, filename):
    """
    Load a dynamic function from a file and attach it to the obj. (A Network
    or Trajectory.)
    
    The filename must be <function_name>.py
    """
    f = file(filename, 'r')
    function_body = f.read()
    f.close()

    basename = os.path.basename(filename)
    func = os.path.splitext(basename)[0]
    setattr(obj, '%s_functionBody' % func, function_body)
    # We try to get the attribute 'namespace' for the object.
    Network_mod._exec_dynamic_func(obj, func, getattr(obj, 'namespace', {}))

def output_dynamic_functions(obj, directory = SloppyCell._TEMP_DIR):
    """
    Output .py files for this objects's dynamic functions into the given
    directory.
    """
    for func in obj._dynamic_funcs:
        body = getattr(obj, '%s_functionBody' % func, None)
        if body is not None:
            f = file(os.path.join(directory, '%s.py' % func), 'w')
            f.write(getattr(obj, '%s_functionBody' % func))
            f.close()
