import os
import sys
try:
    import Ft
    HAVE_FT = True
    RXNNETS_DIR = os.path.dirname(__file__)
except ImportError:
    HAVE_FT = False
import libsbml
import SloppyCell.ReactionNetworks.Network as Network
import SloppyCell.ReactionNetworks.Reactions as Reactions
import SloppyCell.ReactionNetworks.Parsing as Parsing

def toSBMLFile(net, fileName):
    sbmlStr = toSBMLString(net)
    f = file(fileName, 'w')
    f.write(sbmlStr)
    f.close()

def SBMLtoDOT(sbmlFileName, dotFileName):
    if not HAVE_FT:
        print 'Need 4Suite library installed to export to .dot'
        return

    from Ft.Xml.Xslt import Processor
    processor = Processor.Processor()
    from Ft.Xml import InputSource
        
    transformFileName = os.path.join(RXNNETS_DIR, 'sbml_l2v1_todot.xsl')
    transformFile = file(transformFileName, 'r')
    transformURI = 'file://%s' % transformFileName
    transform = InputSource.DefaultFactory.fromStream(transformFile, 
                                                      transformURI)

    processor.appendStylesheet(transform)

    sourceFile = file(sbmlFileName, 'r')
    sourceURI = 'file://%s' % os.path.abspath(sbmlFileName)

    source = InputSource.DefaultFactory.fromStream(sourceFile, sourceURI)

    out = processor.run(source)
    f = file(dotFileName, 'w')
    f.write(out)
    f.close()
    

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
        sc.setConstant(c.isConstant)
        sc.setSize(c.initialValue)
        m.addCompartment(sc)
    
    for id, s in net.species.items():
        ss = libsbml.Species(id)
        ss.setName(s.name)
        ss.setCompartment(s.compartment)
        if s.initialValue is not None:
            ss.setInitialConcentration(s.initialValue)
        ss.setBoundaryCondition(s.isBoundaryCondition)
        m.addSpecies(ss)
    
    for id, p in net.parameters.items():
        sp = libsbml.Parameter(id)
        sp.setName(p.name)
        if p.initialValue is not None:
            sp.setValue(p.initialValue)
        sp.setConstant(p.isConstant)
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

def fromSBMLFile(fileName, id = None):
    f = file(fileName, 'r')
    net = fromSBMLString(f.read(), id)
    f.close()
    return net

def fromSBMLString(sbmlStr, id = None):
    r = libsbml.SBMLReader()
    d = r.readSBMLFromString(sbmlStr)
    m = d.getModel()

    modelId = m.getId()
    if (id == None) and (modelId == ''):
        raise ValueError, 'Network id not specified in SBML or passed in.'
    elif id is not None:
        modelId = id
        
    rn = Network.Network(id = modelId, name = m.getName())

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
                      isBoundaryCondition = isBC,
                      name = name)

    for p in m.getListOfParameters():
        parameter = createNetworkParameter(p)
        rn.addVariable(parameter)

    for rxn in m.getListOfReactions():
        id, name = rxn.getId(), rxn.getName()
        kL = rxn.getKineticLaw()
        kLFormula = kL.getFormula()

        for p in kL.getListOfParameters():
            parameter = createNetworkParameter(p)
            if parameter.id in rn.variables.keys():
                oldId = parameter.id
                parameter.id = id + '_' + parameter.id
                kLFormula = Parsing.\
                        substituteVariableNamesInString(kLFormula, oldId, 
                                                        parameter.id)
            rn.addVariable(parameter)
    
        stoichiometry = {}
        for reactant in rxn.getListOfReactants():
            if reactant.getStoichiometryMath() == None:
                stoichiometry[reactant.getSpecies()] =\
                        -1 * reactant.getStoichiometry()
    
        for product in rxn.getListOfProducts():
            if product.getStoichiometryMath() == None:
                stoichiometry[product.getSpecies()] = product.getStoichiometry()

        for modifier in rxn.getListOfModifiers():
            stoichiometry.setdefault(modifier.getSpecies(), 0)

        rn.addReaction(id = id, stoichiometry = stoichiometry,
                       kineticLaw = kLFormula)

    for ii, r in enumerate(m.getListOfRules()):
        if r.getTypeCode() == libsbml.SBML_ALGEBRAIC_RULE:
            print >> sys.stderr, '*'*20
            print >> sys.stderr, 'Warning: Alegbraic rule specied in SBML file.'
            print >> sys.stderr, 'Algebraic constraints are not implemented in our system!'
            print
            print >> sys.stderr, 'Rule is: %s' % libsbml.formulaToString(r.getMath())
            print >> sys.stderr, '*'*20
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
                    delay = delay, name = name, isTerminal = True)

    return rn

def createNetworkParameter(p):
    id, name = p.getId(), p.getName()
    v = p.getValue()
    isConstant = p.getConstant()

    parameter = Network.Parameter(id = id, value = v, isConstant = isConstant,
                                  name = name)

    return parameter
