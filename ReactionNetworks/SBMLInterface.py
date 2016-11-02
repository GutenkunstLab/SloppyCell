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
import SloppyCell.KeyedList_mod
KeyedList = SloppyCell.KeyedList_mod.KeyedList

# sbml_level and sbml_version are default parameters to pass to
# constructors for libsbml 4.0
sbml_level = 2
sbml_version = 3

def toSBMLFile(net, fileName):
    sbmlStr = toSBMLString(net)
    f = file(fileName, 'w')
    f.write(sbmlStr)
    f.close()

def SBMLtoDOT(sbmlFileName, dotFileName):
    raise DeprecationWarning, 'SBMLtoDOT has been deprecated. Instead, use IO.net_DOT_file(net, filename)'

def rxn_add_stoich(srxn, rid, stoich, is_product=True):
    try:
        stoich = float(stoich)
        if stoich < 0:
            sr = srxn.createReactant()
            sr.setStoichiometry(-stoich)
        elif stoich > 0:
            sr = srxn.createProduct()
            sr.setStoichiometry(stoich)
        elif stoich == 0:
            sr = srxn.createModifier()
        sr.setSpecies(rid)

    except ValueError:
        formula = stoich.replace('**', '^')
        math_ast = libsbml.parseFormula(formula)
        #try:
        #    smath = libsbml.StoichiometryMath(math_ast)
        #except:
        # normally this is in an except block but the try above doesn't throw erro
        # I am expecting in libSBML 4.0
        try:
            smath = libsbml.StoichiometryMath(math_ast)
        except NotImplementedError:
            smath = libsbml.StoichiometryMath(sbml_level, sbml_version)
            smath.setMath(math_ast)

        if is_product == True:
            sr = srxn.createProduct()
        else:
            sr = srxn.createReactant()

        sr.setSpecies(rid)
        sr.setStoichiometryMath(smath)

def replaceTime(ast):
    """
    Takes an ast and recurively checks for any variables with id
    't' or 'time'.  Then, changes the type of those ast nodes with
    AST_NAME_TIME so that they get output in SBML with the proper
    MathML designation.

    Based on code proposed by Darren Wilkinson on the libsbml
    forum (http://bit.ly/8Yl5sx) .
    """
    if (ast.getType()==libsbml.AST_NAME):
        if ((ast.getName()=='t') or (ast.getName()=='time')):
            ast.setType(libsbml.AST_NAME_TIME)
    for node in range(ast.getNumChildren()):
        replaceTime(ast.getChild(node))

def toSBMLString(net):
    metaId = 0
    try:
        m = libsbml.Model(net.id)
    except NotImplementedError:
        m = libsbml.Model(sbml_level, sbml_version)
        m.setId(net.id)
    m.setName(net.name)
    m.setMetaId('SloppyCell_{0:05d}'.format(metaId))
    metaId += 1

    for id, fd in net.functionDefinitions.items():
        try:
            sfd = libsbml.FunctionDefinition(id)
        except:
            sfd = libsbml.FunctionDefinition(sbml_level, sbml_version)
            sfd.setId(id)
        sfd.setName(fd.name)
        formula = fd.math
        formula = formula.replace('**', '^')
        formula = 'lambda(%s, %s)' % (','.join(fd.variables), formula)
        sfd.setMath(libsbml.parseFormula(formula))
        sfd.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addFunctionDefinition(sfd)

    for id, c in net.compartments.items():
        try:
            sc = libsbml.Compartment(id)
        except NotImplementedError:
            sc = libsbml.Compartment(sbml_level, sbml_version)
            sc.setId(id)
        sc.setName(c.name)
        sc.setConstant(c.is_constant)
        sc.setSize(c.initialValue)
        sc.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addCompartment(sc)

    for id, s in net.species.items():
        try:
            ss = libsbml.Species(id)
        except NotImplementedError:
            ss = libsbml.Species(sbml_level, sbml_version)
            ss.setId(id)
        ss.setName(s.name)
        ss.setCompartment(s.compartment)
        if s.initialValue is not None and not isinstance(s.initialValue, str):
            ss.setInitialConcentration(s.initialValue)
        ss.setBoundaryCondition(s.is_boundary_condition)
        ss.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addSpecies(ss)

    for id, p in net.parameters.items():
        try:
            sp = libsbml.Parameter(id)
        except NotImplementedError:
            sp = libsbml.Parameter(sbml_level, sbml_version)
            sp.setId(id)
        sp.setName(p.name)
        if p.initialValue is not None:
            sp.setValue(p.initialValue)
        sp.setConstant(p.is_constant)
        sp.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addParameter(sp)

    for id, r in net.rateRules.items():
        try:
            sr = libsbml.RateRule()
        except NotImplementedError:
            sr = libsbml.RateRule(sbml_level, sbml_version)
        sr.setVariable(id)
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        sr.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addRule(sr)

    for id, r in net.assignmentRules.items():
        try:
            sr = libsbml.AssignmentRule()
        except NotImplementedError:
            sr = libsbml.AssignmentRule(sbml_level, sbml_version)
        sr.setVariable(id)
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        sr.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addRule(sr)

    for r, r in net.algebraicRules.items():
        try:
            sr = libsbml.AlgebraicRule()
        except NotImplementedError:
            sr = libsbml.AlgebraicRule(sbml_level, sbml_version)
        formula = r.replace('**', '^')
        sr.setMath(libsbml.parseFormula(formula))
        sr.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addRule(sr)

    for id, rxn in net.reactions.items():
        # Need to identify modifiers in kinetic law and add them to
        # stoichiometry
        kl_vars = ExprManip.extract_vars(rxn.kineticLaw)
        species_in_kl = kl_vars.intersection(net.species.keys())
        for s in species_in_kl:
            if not rxn.stoichiometry.has_key(s):
                rxn.stoichiometry[s] = 0

        try:
            srxn = libsbml.Reaction(id)
        except NotImplementedError:
            srxn = libsbml.Reaction(sbml_level, sbml_version)
            srxn.setId(id)
        srxn.setName(rxn.name)
        # Handle the case where the model was originally read in from an
        # SBML file, so that the reactants and products of the Reaction
        # object are explicitly set.
        if rxn.reactant_stoichiometry != None and \
            rxn.product_stoichiometry != None:
            for rid, stoich_list in rxn.reactant_stoichiometry.items():
                for stoich in stoich_list:
                    rxn_add_stoich(srxn, rid, -float(stoich), is_product=False)
            for rid, stoich_list in rxn.product_stoichiometry.items():
                for stoich in stoich_list:
                    rxn_add_stoich(srxn, rid, stoich, is_product=True)
        # Handle the case where the model was created using the SloppyCell
        # API, in which case reactants and products are inferred from their
        # stoichiometries
        else:
            for rid, stoich in rxn.stoichiometry.items():
                rxn_add_stoich(srxn, rid, stoich)

        formula = rxn.kineticLaw.replace('**', '^')
        try:
            kl = libsbml.KineticLaw(formula)
        except NotImplementedError:
            kl = libsbml.KineticLaw(sbml_level, sbml_version)
            kl.setFormula(formula)
        srxn.setKineticLaw(kl)
        srxn.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addReaction(srxn)

    for id, e in net.events.items():
        try:
            se = libsbml.Event(id)
        except NotImplementedError:
            se = libsbml.Event(sbml_level, sbml_version)
            se.setId(id)
        se.setName(e.name)
        formula = e.trigger.replace('**', '^')
        formula = formula.replace('and_func(', 'and(')
        formula = formula.replace('or_func(', 'or(')

        ast = libsbml.parseFormula(formula)
        if ast is None:
            raise ValueError('Problem parsing event trigger: %s. Problem may '
                             'be use of relational operators (< and >) rather '
                             'than libsbml-friendly functions lt and gt.'
                             % formula)
        try:
            se.setTrigger(ast)
        except TypeError:
            try:
                trigger = libsbml.Trigger(ast)
            except NotImplementedError:
                trigger = libsbml.Trigger(sbml_level, sbml_version)
                trigger.setMath(ast)
            se.setTrigger(trigger)
        formula = str(e.delay).replace('**', '^')
        try:
            se.setDelay(libsbml.parseFormula(formula))
        except TypeError:
            try:
                se.setDelay(libsbml.Delay(libsbml.parseFormula(formula)))
            except NotImplementedError:
                delay = libsbml.Delay(sbml_level, sbml_version)
                delay.setMath(libsbml.parseFormula(formula))
                se.setDelay(delay)
        for varId, formula in e.event_assignments.items():
            try:
                sea = libsbml.EventAssignment()
            except NotImplementedError:
                sea = libsbml.EventAssignment(sbml_level, sbml_version)
            sea.setVariable(varId)
            formula = str(formula).replace('**', '^')
            formula = formula.replace('and_func(', 'and(')
            formula = formula.replace('or_func(', 'or(')
            ast = libsbml.parseFormula(formula)
            replaceTime(ast)
            sea.setMath(ast)
            se.addEventAssignment(sea)
        se.setMetaId('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1
        m.addEvent(se)

    for id, con in net.constraints.items():
        try:
            scon = libsbml.Constraint()
        except NotImplementedError:
            scon = libsbml.Constraint(sbml_level, sbml_version)
        scon.setId(con.id)
        scon.setName(con.name)
        formula = con.trigger.replace('**', '^')
        ast = libsbml.parseFormula(formula)
        if ast is None:
            raise ValueError('Problem parsing constraint math: %s. Problem may '
                             'be use of relational operators (< and >) rather '
                             'than libsbml-friendly functions lt and gt.'
                             % formula)
        scon.setMath(ast)
        se.setcon('SloppyCell_{0:05d}'.format(metaId))
        metaId += 1

        m.addConstraint(scon)

    d = libsbml.SBMLDocument(sbml_level, sbml_version)
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

def fromSBMLFile(fileName, id = None, duplicate_rxn_params=False):
    f = file(fileName, 'r')
    net = fromSBMLString(f.read(), id, duplicate_rxn_params)
    f.close()
    return net

def stoichToString(species, stoich):
    if stoich is None:
        stoich = str(species.getStoichiometry())
    elif hasattr(stoich, 'getMath'): # libsbml > 3.0
        stoich = libsbml.formulaToString(stoich.getMath())
    else: # libsbml 2.3.4
        stoich = libsbml.formulaToString(stoich)
    return stoich

def fromSBMLString(sbmlStr, id = None, duplicate_rxn_params=False):
    r = libsbml.SBMLReader()
    d = r.readSBMLFromString(sbmlStr)
    if d.getNumErrors():
        message = 'libSBML reported errors in SBML file. Try running file '\
                'through the online validator: '\
                'http://www.sbml.org/Facilities/Validator . Specific errors '\
                'noted are: '
        errors = []
        for ii in range(d.getNumErrors()):
            pm = d.getError(ii)
            errors.append(pm.getMessage())
        raise ValueError(message + '; '.join(errors))

    m = d.getModel()

    modelId = m.getId()
    if (id is None) and (modelId == ''):
        raise ValueError('Network id not specified in SBML or passed in.')
    elif id is not None:
        modelId = id

    rn = Network_mod.Network(id = modelId, name = m.getName())

    for f in m.getListOfFunctionDefinitions():
        id, name = f.getId(), f.getName()
        math = f.getMath()
        variables = []
        for ii in range(math.getNumChildren() - 1):
            variables.append(libsbml.formulaToString(math.getChild(ii)))

        math = libsbml.formulaToString(math.getRightChild())

        rn.addFunctionDefinition(id, variables, math)

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
        if s.isSetInitialConcentration():
            iC = s.getInitialConcentration()
        elif s.isSetInitialAmount():
            iC = s.getInitialAmount()
        else:
            iC = 1
        isBC, isConstant = s.getBoundaryCondition(), s.getConstant()

        xml_text = s.toSBML()
        uniprot_ids = set([entry[1:].split('"')[0]
                           for entry in xml_text.split('uniprot')[1:]])

	rn.addSpecies(id = id, compartment = compartment,
                      initialConcentration = iC,
                      isConstant = isConstant,
                      is_boundary_condition = isBC,
                      name = name, uniprot_ids = uniprot_ids)

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
        # For each species appearing in the reaction, we build up a string
        # representing the stoichiometry. Then we'll simplify that string and
        # see whether we ended up with a float value in the end.
        stoichiometry = {}
        reactant_stoichiometry = {}
        product_stoichiometry = {}
        for reactant in rxn.getListOfReactants():
            species = reactant.getSpecies()
            stoichiometry.setdefault(species, '0')
            stoich = reactant.getStoichiometryMath()
            stoich = stoichToString(reactant, stoich)
            stoichiometry[species] += '-(%s)' % stoich
            if species in reactant_stoichiometry:
                reactant_stoichiometry[species].append(stoich)
            else:
                reactant_stoichiometry[species] = [stoich]

        for product in rxn.getListOfProducts():
            species = product.getSpecies()
            stoichiometry.setdefault(species, '0')
            stoich = product.getStoichiometryMath()
            stoich = stoichToString(product, stoich)
            stoichiometry[species] += '+(%s)' % stoich
            if species in product_stoichiometry:
                product_stoichiometry[species].append(stoich)
            else:
                product_stoichiometry[species] = [stoich]

        for species, stoich in stoichiometry.items():
            stoich = ExprManip.simplify_expr(stoich)
            try:
                # Try converting the string to a float.
                stoich = float(stoich)
            except ValueError:
                pass
            stoichiometry[species] = stoich

        for modifier in rxn.getListOfModifiers():
            stoichiometry.setdefault(modifier.getSpecies(), 0)

        rn.addReaction(id = id, stoichiometry = stoichiometry,
                       kineticLaw = kLFormula,
                       reactant_stoichiometry = reactant_stoichiometry,
                       product_stoichiometry = product_stoichiometry)

    for ii, r in enumerate(m.getListOfRules()):
        if r.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE:
            math = libsbml.formulaToString(r.getMath())
            rn.add_algebraic_rule(math)
        else:
            variable = r.getVariable()
            math = libsbml.formulaToString(r.getMath())
            if r.getTypeCode() == libsbml.SBML_ASSIGNMENT_RULE:
                rn.addAssignmentRule(variable, math)
            elif r.getTypeCode() == libsbml.SBML_RATE_RULE:
                rn.addRateRule(variable, math)



    for ii, e in enumerate(m.getListOfEvents()):
        id, name = e.getId(), e.getName()

        if id == '':
            id = 'event%i' % ii

        try:
            # For libSBML 3.0
            trigger_math = e.getTrigger().getMath()
        except AttributeError:
            # For older versions
            trigger_math = e.getTrigger()
        trigger = libsbml.formulaToString(trigger_math)
        trigger = trigger.replace('or(','or_func(')
        trigger = trigger.replace('and(','and_func(')

        if e.getDelay() is not None:
            try:
                # For libSBML 3.0
                delay_math = e.getDelay().getMath()
            except AttributeError:
                # For older versions
                delay_math = e.getDelay()
            delay = libsbml.formulaToString(delay_math)
        else:
            delay = 0

        timeUnits = e.getTimeUnits()
        eaDict = KeyedList()
        for ea in e.getListOfEventAssignments():
            ea_formula = libsbml.formulaToString(ea.getMath())
            ea_formula = ea_formula.replace('or(','or_func(')
            ea_formula = ea_formula.replace('and(','and_func(')
            eaDict.set(ea.getVariable(), ea_formula)

        rn.addEvent(id = id, trigger = trigger, eventAssignments = eaDict,
                    delay = delay, name = name)


    for ii, con in enumerate(m.getListOfConstraints()):
        id, name = con.getId(), con.getName()
        if id == '':
            id = 'constraint%i' % ii

        trigger_math = con.getMath()

        trigger = libsbml.formulaToString(trigger_math)

        if con.isSetMessage():
            message = con.getMessage()
        else:
            message = None

        rn.addConstraint(id = id, trigger = trigger, message = message,
                    name = name)

    return rn

def createNetworkParameter(p):
    id, name = p.getId(), p.getName()
    v = p.getValue()
    isConstant = p.getConstant()

    parameter = Network_mod.Parameter(id = id, value = v, is_constant = isConstant,
                                      name = name, typical_value = None,
                                      is_optimizable = True)
                                      # optimizable by default

    return parameter

