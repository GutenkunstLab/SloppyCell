import os
import sets

try:
    import SBMLInterface as SBML
    to_SBML_file = SBML.toSBMLFile
    from_SBML_file = SBML.fromSBMLFile
except ImportError:
    print 'Failed to import libsbml.'
    print 'SBML import and export not available.'

import SloppyCell
from SloppyCell.KeyedList_mod import KeyedList
import SloppyCell.ExprManip as ExprManip
import Network_mod
import Trajectory_mod

def net_DOT_file(net, filename = None, size=(7.5,10)):
    lines = []
    lines.append('digraph "%s" {' % net.id)
    lines.append('\tsize="'+str(size)+'!"')
    lines.append('\tratio=fill')
    for id in net.species.keys():
        lines.append('\t"%s"[color=black]' % net.get_component_name(id))

    lines.append('')
    for reac_id, rxn in net.reactions.items():
        rxn_name = net.get_component_name(reac_id)
        lines.append('\t"%s"[shape=box][color=red]' % rxn_name)
        for rid, stoich in rxn.stoichiometry.items():
            rname = net.get_component_name(rid)
            if stoich > 0:
                lines.append('\t\t"%s" -> "%s";' % (rxn_name, rname))
            elif stoich < 0:
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


def eqns_TeX_file(net, filename=None, simpleTeX=False, landscape=False):
    """
    Output a TeX file containing the network equations and other information.

    net: Network to work with.
    filename: Filename for resulting TeX file. If filename==None, the output 
              file will be <network id>.tex
    simpleTeX: If True, some TeX that causes problems in WYSIWYG editors such as
               LyX will be omitted.
    """
    net.compile()

    lines = []
    lines.append(r'\documentclass{article}')
    lines.append(r'\usepackage{amsmath}')
    lines.append(r'\usepackage{fullpage}')
    lines.append(r'\usepackage{longtable}')
    if landscape == True:
        lines.append(r'\usepackage[a4paper,landscape]{geometry}')
    lines.append(r'\begin{document}')
    lines.append(_net_eqns_to_TeX(net, simpleTeX))
    lines.append(r'\end{document}')

    if filename is None:
        filename = '%s.tex' % net.id
    f = file(filename, 'w')
    f.write(os.linesep.join(lines))
    f.close()

def _net_eqns_to_TeX(net, simpleTeX=False):
    """
    Return a string that contains the longtable-bound TeX'd equations for the
    network. Also includes events and current optimizable parameter values.

    simpleTeX: If True, some TeX that causes problems in WYSIWYG editors such as
               LyX will be omitted.
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

    outputs = []
    if net.functionDefinitions:
        func_KL = KeyedList()
        for func_id, func in net.functionDefinitions.items():
            lhs = '%s(%s)' % (func_id, r','.join(func.variables))
            func_KL.set(lhs, func.math)
        funcs_str = ExprManip.Py2TeX.dict2TeX(func_KL, name_dict, 
                                              split_terms=False,
                                              simpleTeX=simpleTeX)
        outputs.extend([r'\section*{Function Definitions}', funcs_str])

    if net.assignmentRules:
        assigns_str = ExprManip.Py2TeX.dict2TeX(net.assignmentRules, name_dict, 
                                                split_terms=True,
                                                simpleTeX=simpleTeX)
        outputs.extend([r'\section*{Assignment Rules}', assigns_str])

    if net.diff_eq_rhs:
        eqns_str = ExprManip.Py2TeX.dict2TeX(net.diff_eq_rhs, name_dict, 
                                             r'\frac{d\,%s}{dt}', 
                                             split_terms=True,
                                             simpleTeX=simpleTeX)
        outputs.extend([r'\section*{Differential Equations}', eqns_str])

    if net.algebraicRules:
        eqns_str = ExprManip.Py2TeX.dict2TeX(net.algebraicRules, name_dict, 
                                             lhs_form='0', split_terms=True,
                                             simpleTeX=simpleTeX)
        outputs.extend([r'\section*{Algebraic Equations}', eqns_str])

    if net.events:
        outputs.append(r'\section*{Events}')

        for e in net.events:
            e_name = net.get_component_name(e.id, TeX_form=True)
            outputs.append(r'$%s$' % (e_name))
            trigger_str = ExprManip.Py2TeX.expr2TeX(e.trigger, name_dict)
            outputs.append(r'Trigger: $%s$' % trigger_str)
            
            for id, result in e.event_assignments.items():
                outputs.append(r'\begin{itemize}')
                id_str = ExprManip.Py2TeX.expr2TeX(id, name_dict)
                rule_str = ExprManip.Py2TeX.expr2TeX(str(result), name_dict)
                outputs.append(r'\item $%s = %s$' % (id_str, rule_str))
                outputs.append(r'\end{itemize}')





    if net.optimizableVars:
        outputs.append(r'\section*{Optimizable Parameters}')
        params = net.GetParameters()
        outputs.append(r'\begin{longtable}{|r|l|}')
        outputs.append(r'\hline')
        for id, value in params.items():
            id_str = ExprManip.Py2TeX.expr2TeX(id, name_dict)
            value_str = ExprManip.Py2TeX.expr2TeX(str(value), name_dict)
            outputs.append(r'$%s$& $%s$\\' % (id_str, value_str))
            outputs.append(r'\hline')
        outputs.append(r'\end{longtable}')

    outputs_copy = []
    for line in outputs:
        # fix and_func and or_func
        line = line.replace('and_func','and\_func')
        line = line.replace('or_func','or\_func')
        outputs_copy.append(line)

    return os.linesep.join(outputs_copy)
                            

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
    file_type = os.path.splitext(basename)[1][1:]

    if isinstance(obj, Trajectory_mod.Trajectory):
        setattr(obj, '%s_functionBody' % func, function_body)
        # We try to get the attribute 'namespace' for the object.
        Network_mod._exec_dynamic_func(obj, func, getattr(obj, 'namespace', {}))
    elif isinstance(obj, Network_mod.Network):
        if file_type == 'py':
            obj._dynamic_funcs_python[func] = function_body
            obj.exec_dynamic_functions(disable_c = True)
        elif file_type == 'c':
            obj.exec_dynamic_functions(curr_c_code = function_body)

def output_dynamic_functions(obj, directory = SloppyCell._TEMP_DIR):
    """
    Output .py files for this objects's dynamic functions into the given
    directory.
    """
    if isinstance(obj, Trajectory_mod.Trajectory):
        for func in obj._dynamic_funcs:
            body = getattr(obj, '%s_functionBody' % func, None)
            if body is not None:
                f = file(os.path.join(directory, '%s.py' % func), 'w')
                f.write(getattr(obj, '%s_functionBody' % func))
                f.close()
    elif isinstance(obj, Network_mod.Network):
        for func, body in obj._dynamic_funcs_python.items():
            if body != None:
                f = file(os.path.join(directory, '%s.py' % func), 'w')
                f.write(body)
                f.close()
        c_code = obj.get_c_code()
        f = file(os.path.join(directory, '%s.c' % obj.get_id()), 'w')
        f.write(c_code)
        f.close()
