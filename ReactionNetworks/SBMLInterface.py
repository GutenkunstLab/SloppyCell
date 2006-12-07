"""
Methods for loading from and saving to SBML files.
"""
__docformat__ = "restructuredtext en"

import os
import sys
import logging
logger = logging.getLogger('RxnNets.SBMLInterface')

import libsbml

import Network_mod
import SloppyCell.ExprManip as ExprManip

def toSBMLFile(net, fileName):
    sbmlStr = toSBMLString(net)
    f = file(fileName, 'w')
    f.write(sbmlStr)
    f.close()

def SBMLtoDOT(sbmlFileName, dotFileName):
    raise DeprecationWarning, 'SBMLtoDOT has been deprecated. Instead, use IO.net_DOT_file(net, filename)'

def toSBMLString(net):
    m = libsbml.Model(net.id)
    m.setName(net.name)
    
    for id, fd in net.functionDefinitions.items():
        sfd = libsbml.FunctionDefinition(id)
        sfd.setName(fd.name)
        formula = fd.math
        formula = formula.replace('**', '^')
        formula = 'lambda(%s, %s)' % (','.join(fd.variables), formula)
        sfd.setMath(libsbml.parseFormula(formula))
        m.addFunctionDefinition(sfd)
    
    for id, c in net.compartments.items():
        sc = libsbml.Compartment(id)
        sc.setName(c.name)
        sc.setConstant(c.is_constant)
        sc.setSize(c.initialValue)
        m.addCompartment(sc)
    
    for id, s in net.species.items():
        ss = libsbml.Species(id)
        ss.setName(s.name)
        ss.setCompartment(s.compartment)
        if s.initialValue is not None and not isinstance(s.initialValue, str):
            ss.setInitialConcentration(s.initialValue)
        ss.setBoundaryCondition(s.is_boundary_condition)
        m.addSpecies(ss)
    
    for id, p in net.parameters.items():
        sp = libsbml.Parameter(id)
        sp.setName(p.name)
        if p.initialValue is not None:
            sp.setValue(p.initialValue)
        sp.setConstant(p.is_constant)
        m.addParameter(sp)

    for id, r in net.rateRules.items():
        sr = libsbml.RateRule()
        sr.setVariable(id)
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        m.addRule(sr)

    for id, r in net.assignmentRules.items():
        sr = libsbml.AssignmentRule()
        sr.setVariable(id)
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        m.addRule(sr)

    for r, r in net.algebraicRules.items():
        sr = libsbml.AlgebraicRule()
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        m.addRule(sr)
        
    for id, rxn in net.reactions.items():
        srxn = libsbml.Reaction(id)
        srxn.setName(rxn.name)
        for rid, stoich in rxn.stoichiometry.items():
            sr = libsbml.SpeciesReference(rid)
            if isinstance(stoich, str):
                formula = stoich.replace('**', '^')
                sr.setStoichiometryMath(libsbml.parseFormula(formula))
                srxn.addReactant(sr)
            else:
                if stoich < 0:
                    sr.setStoichiometry(-stoich)
                    srxn.addReactant(sr)
                if stoich > 0:
                    sr.setStoichiometry(stoich)
                    srxn.addProduct(sr)
                if stoich == 0:
                    sr = libsbml.ModifierSpeciesReference(rid)
                    srxn.addModifier(sr)
        formula = rxn.kineticLaw.replace('**', '^')
        kl = libsbml.KineticLaw(formula)
        srxn.setKineticLaw(kl)
        m.addReaction(srxn)
    
    for id, e in net.events.items():
        se = libsbml.Event()
        se.setName(e.name)
        formula = e.trigger.replace('**', '^')
        se.setTrigger(libsbml.parseFormula(formula))
        formula = str(e.delay).replace('**', '^')
        se.setDelay(libsbml.parseFormula(formula))
        for varId, formula in e.eventAssignments.items():
            sea = libsbml.EventAssignment()
            sea.setVariable(varId)
            formula = str(formula).replace('**', '^')
            sea.setMath(libsbml.parseFormula(formula))
            se.addEventAssignment(sea)
        m.addEvent(se)
    
    d = libsbml.SBMLDocument()
    d.setModel(m)
    sbmlStr = libsbml.writeSBMLToString(d)

    return sbmlStr

def _print_sbml_fatals(doc):
    for ii in range(doc.getNumFatals()):
        logger.critical(d.getFatal(ii).getMessage())

def to_SBML_l2v1(from_name, to_name):
    """
    Convert an SBML file to level 2, version 1 using lisbml.

    from_name  Name of file to read SBML from.
    to_name    File to output SBML to.
    """
    doc = libsbml.readSBML(os.path.abspath(from_name))
    if isinstance(doc, int):
        logger.critical('Fatal Errors reading SBML from file %s!' % from_name)
        _print_sbml_fatals(doc)
    errors = doc.setLevel(2)
    if errors:
        logger.critical('Fatal Errors converting %f to level 2!' % from_name)
        _print_sbml_fatals(doc)
    errors = doc.setVersion(1)
    if errors:
        logger.critical('Fatal Errors converting %f to level 2, version 1!'
                        % from_name)
        _print_sbml_fatals(doc)
    success = libsbml.writeSBML(doc, to_name)
    if not success:
        logger.critical('Error writing to %s' % to_name)

def fromSBMLFile(fileName, id = None):
    f = file(fileName, 'r')
    net = fromSBMLString(f.read(), id)
    f.close()
    return net

def fromSBMLString(sbmlStr, id = None, duplicate_rxn_params=False):
    r = libsbml.SBMLReader()
    d = r.readSBMLFromString(sbmlStr)
    m = d.getModel()

    modelId = m.getId()
    if (id == None) and (modelId == ''):
        raise ValueError, 'Network id not specified in SBML or passed in.'
    elif id is not None:
        modelId = id
        
    rn = Network_mod.Network(id = modelId, name = m.getName())

    for c in m.getListOfCompartments():
        id, name = c.getId(), c.getName()
        size = c.getSize()
        isConstant = c.getConstant()

        rn.addCompartment(id = id, size = size, 
                          isConstant = isConstant, 
                          name = name)

    for s in m.getListOfSpecies():
        id, name = s.getId(), s.getName()
        compartment = s.getCompartment()
        iC = s.getInitialConcentration()
        isBC, isConstant = s.getBoundaryCondition(), s.getConstant()
	
	rn.addSpecies(id = id, compartment = compartment,
                      initialConcentration = iC,
                      isConstant = isConstant,
                      is_boundary_condition = isBC,
                      name = name)
        #rn.addSpecies(id = id, compartment = compartment,
        #              initial_conc = iC,
        #              is_constant = isConstant,
        #              is_boundary_condition = isBC,
        #              name = name)

    for p in m.getListOfParameters():
        parameter = createNetworkParameter(p)
        rn.addVariable(parameter)

    for rxn in m.getListOfReactions():
        id, name = rxn.getId(), rxn.getName()
        kL = rxn.getKineticLaw()
        kLFormula = kL.getFormula()

        substitution_dict = {}
        # Deal with parameters defined within reactions
        for p in kL.getListOfParameters():
            parameter = createNetworkParameter(p)
            # If a parameter with this name already exists, **and it has a
            # different value than this parameter** we rename this parameter
            # instance by prefixing it with the rxn name so there isn't a
            # clash.
            if parameter.id in rn.variables.keys():
                logger.warn('Parameter %s appears in two different reactions '
                            'in SBML file.' % parameter.id)
                if parameter.value != rn.variables.get(parameter.id).value or\
                   duplicate_rxn_params:
                    oldId = parameter.id
                    parameter.id = id + '_' + parameter.id
                    substitution_dict[oldId] = parameter.id
                    logger.warn('It has different values in the two positions '
                                'so we are creating a new parameter %s.'
                                % (parameter.id))
                else:
                    logger.warn('It has the same value in the two positions '
                                'so we are only defining one parameter %s. '
                                'This behavior can be changed with the option '
                                'duplicate_rxn_params = True' % (parameter.id))

            if parameter.id not in rn.variables.keys():
                rn.addVariable(parameter)
        kLFormula = ExprManip.sub_for_vars(kLFormula, substitution_dict) 
    
        # Assemble the stoichiometry. SBML has the annoying trait that 
        #  species can appear as both products and reactants and 'cancel out'
        stoichiometry = {}
        for reactant in rxn.getListOfReactants():
            stoichiometry.setdefault(reactant.getSpecies(), 0)
            if reactant.getStoichiometryMath() == None:
                stoichiometry[reactant.getSpecies()] -=\
                        reactant.getStoichiometry()
    
        for product in rxn.getListOfProducts():
            stoichiometry.setdefault(product.getSpecies(), 0)
            if product.getStoichiometryMath() == None:
                stoichiometry[product.getSpecies()] += product.getStoichiometry()

        for modifier in rxn.getListOfModifiers():
            stoichiometry.setdefault(modifier.getSpecies(), 0)

        rn.addReaction(id = id, stoichiometry = stoichiometry,
                       kineticLaw = kLFormula)

    for ii, r in enumerate(m.getListOfRules()):
        if r.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE:
            print >> sys.stderr, '*'*20
            print >> sys.stderr, 'Warning: Alegbraic rule specied in SBML file.'
            print >> sys.stderr, 'Algebraic constraints are still being tested in our system!'
            print
            print >> sys.stderr, 'Rule is: %s' % libsbml.formulaToString(r.getMath())
            print >> sys.stderr, '*'*20
            math = libsbml.formulaToString(r.getMath())
            rn.add_algebraic_rule(math)
        else:
            variable = r.getVariable()
            math = libsbml.formulaToString(r.getMath())
            if r.getTypeCode() == libsbml.SBML_ASSIGNMENT_RULE:
                rn.addAssignmentRule(variable, math)
            elif r.getTypeCode() == libsbml.SBML_RATE_RULE:
                rn.addRateRule(variable, math)

    for f in m.getListOfFunctionDefinitions():
        id, name = f.getId(), f.getName()
        math = f.getMath()
        variables = []
        for ii in range(math.getNumChildren() - 1):
            variables.append(libsbml.formulaToString(math.getChild(ii)))

        math = libsbml.formulaToString(math.getRightChild())

        rn.addFunctionDefinition(id, variables, math)

    for ii, e in enumerate(m.getListOfEvents()):
        id, name = e.getId(), e.getName()
        if id == '':
            id = 'event%i' % ii

        trigger = libsbml.formulaToString(e.getTrigger())
        if e.getDelay() is not None:
            delay = libsbml.formulaToString(e.getDelay())
        else:
            delay = 0

        timeUnits = e.getTimeUnits()
        eaDict = {}
        for ea in e.getListOfEventAssignments():
            eaDict[ea.getVariable()] = libsbml.formulaToString(ea.getMath())

        rn.addEvent(id = id, trigger = trigger, eventAssignments = eaDict, 
                    delay = delay, name = name)

    return rn

def createNetworkParameter(p):
    id, name = p.getId(), p.getName()
    v = p.getValue()
    isConstant = p.getConstant()

    parameter = Network_mod.Parameter(id = id, value = v, is_constant = isConstant,
                                      name = name, typical_value = None, is_optimizable = True)
				  # optimizable by default

    return parameter

