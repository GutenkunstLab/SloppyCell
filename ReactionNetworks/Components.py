import logging
logger = logging.getLogger('ReactionNetworks.Components')

import sets

import SloppyCell.ExprManip as ExprManip

#
# Containers for the more complex SBML entities.
#

class FunctionDefinition:
    def __init__(self, id, variables, math, name = ''):
        self.id = id
        self.variables = variables
        self.math = math
        self.name = name

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
                (self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not (self == other)
    
class Variable:
    def __init__(self, id, value, 
                 name, typicalValue, 
                 is_constant, is_optimizable):

        if typicalValue is None:
            if value and not isinstance(value, str):
                typicalValue = abs(value)
            else:
                typicalValue = 1

        self.id, self.name = id, name
        self.value, self.initialValue = value, value
        self.typicalValue = typicalValue
        self.is_constant, self.is_optimizable = is_constant, is_optimizable

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
                (self.__dict__ == other.__dict__)

class Compartment(Variable):
    def __init__(self, id, initial_size, name, typical_value, 
                 is_constant, is_optimizable):
        Variable.__init__(self, id, initial_size, name, typical_value, 
                          is_constant, is_optimizable)

class Species(Variable):
    def __init__(self, id, compartment, initial_conc, 
                 name, typical_value,
                 is_boundary_condition, 
                 is_constant, is_optimizable, uniprot_ids=None):
        self.compartment = compartment
        self.is_boundary_condition = is_boundary_condition
        self.uniprot_ids = uniprot_ids
        if uniprot_ids is None:
            self.uniprot_ids = set()

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
    def __init__(self, id, trigger, event_assignments, delay, name,
                 buffer):
        self.id, self.name = id, name
        self.delay = delay 
        self.is_terminal = (len(event_assignments) > 0)

        self.timeTriggered = False
        self.parseTrigger(trigger)

        self.event_assignments = event_assignments
        self.new_values = {}

        self.buffer=buffer
        if (self.buffer > 0 and self.delay != 0):
            logger.warn('Event %s has buffer > 0 and delay != 0. This case '
                        'has not been tested.')
                        

    def __eq__(self, other):
        # This is a little tricky, because Event objects get modified when their
        #  corresponding event fires. (Probably a poor design.)
        # So we need to only compare the important attributes.
        attrs_to_compare = ['__class__', 'id', 'trigger', 'event_assignments',
                            'delay', 'name']
        for attr in attrs_to_compare:
            if getattr(self, attr) != getattr(other, attr):
                return False
        return True

    def __ne__(self, other):
        return not (self == other)
        
    def parseTrigger(self, trigger):
        if '<' in trigger or '>' in trigger or '=' in trigger:
            raise ValueError('Event triggers must use the functions gt and lt, '
                             'rather than the symbols > and <. For example, '
                             'to trigger when B becomes less than A, use '
                             'lt(B,A).')

        # Figures out if the event is time-triggered and parses it to niceness.

        # 'and' is a reserved keyword in python, so the parser will break
        #  unless we substitute the name here.
        trigger = trigger.replace('and(', 'and_func(')

        self.trigger = trigger

        if ExprManip.extract_vars(trigger) == sets.Set(['time']):
            self.timeTriggered = True
            ast = ExprManip.AST.strip_parse(trigger)
            firstArg = ExprManip.AST.ast2str(ast.args[0])
            secondArg = ExprManip.AST.ast2str(ast.args[1])

            if firstArg == 'time':
                self.triggeringTime = eval(secondArg)
            elif secondArg == 'time':
                self.triggeringTime = eval(firstArg)
            else:
                raise 'Problem in time triggered events'

class ConstraintEvent(Event):
    def __init__(self, id, trigger, message, name):
        self.id, self.name = id, name
        self.delay = 0.0
        self.is_terminal = True
        self.event_assignments = {}
        self.message = message

        self.timeTriggered = False
        self.parseTrigger(trigger)
        self.buffer=0.0





# This is a placeholder object that will be used to store information about
#  events during integrations.
class event_info(object):
    pass
