import sets

import symbolic
import Parsing

#
# Containers for the more complex SBML entities.
#

class FunctionDefinition:
    def __init__(self, id, variables, math, name = ''):
        self.id = id
        self.variables = variables
        self.math = math
        self.name = name
    
class Variable:
    def __init__(self, id, value, 
                 name, typicalValue, 
                 isConstant, isOptimizable):
        if typicalValue is None:
            if value != 0:
                typicalValue = value
            else:
                typicalValue = 1

        self.id, self.name = id, name
        self.value, self.initialValue = value, value
        self.typicalValue = typicalValue
        self.isConstant, self.isOptimizable = isConstant, isOptimizable

class Compartment(Variable):
    def __init__(self, id, size, name, is_constant, 
                 typical_value, is_optimizable):
        Variable.__init__(self, id, size, name, typical_value, 
                          is_constant, is_optimizable)

class Species(Variable):
    def __init__(self, id, compartment, initial_conc, 
                 name, typical_value,
                 is_boundary_condition, 
                 is_constant, is_optimizable):
        self.compartment = compartment
        self.is_boundary_condition = is_boundary_condition

        Variable.__init__(self, id, initial_conc, 
                          name, typical_value,
                          is_constant, is_optimizable)

class Parameter(Variable):
    def __init__(self, id, value, name,
                 is_constant, typical_value, is_optimizable):
        Variable.__init__(self, id, value, 
                          name, typical_value,
                          is_constant, is_optimizable)

class Event:
    def __init__(self, id, trigger, eventAssignments, delay = 0, name = '', 
                 isTerminal = True):
        self.id, self.name = id, name
        self.delay = delay 
        self.isTerminal = isTerminal

        self.timeTriggered = False
        self.parseTrigger(trigger)

        self.eventAssignments = eventAssignments
        self.newValues = {}
        
    def parseTrigger(self, trigger):
        # Figures out if the event is time-triggered and parses it to niceness.
        self.trigger = trigger
        ast = symbolic.string2ast(trigger)
        if Parsing.extractVariablesFromAST(ast) == sets.Set(['time']):
            self.timeTriggered = True
            firstArg = symbolic.ast2string(ast[2][2][1])
            secondArg = symbolic.ast2string(ast[2][2][3])

            if firstArg == 'time':
                self.triggeringTime = eval(secondArg)
            elif secondArg == 'time':
                self.triggeringTime = eval(firstArg)
            else:
                raise 'Problem in time triggered events'
