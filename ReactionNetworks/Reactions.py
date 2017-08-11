import sets
import SloppyCell.ExprManip as ExprManip

class Reaction:
    def __init__(self, id, stoichiometry, kineticLaw = '', name = '',
                 reactant_stoichiometry = None, product_stoichiometry = None):
        self.id = id
        self.stoichiometry = stoichiometry
        # self.reactant_stoichiometry and self.product_stoichiometry
        # are defined to help repserve the stoichiometry defined in an
        # SBML model
        self.reactant_stoichiometry = reactant_stoichiometry
        self.product_stoichiometry = product_stoichiometry
        self.kineticLaw = kineticLaw
        self.name = name

        variables = ExprManip.extract_vars(kineticLaw)
        self.parameters = variables.difference(sets.Set(stoichiometry.keys()))

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
                self.kineticLaw == other.kineticLaw and \
                self.stoichiometry == other.stoichiometry

    def __ne__(self, other):
        return not (self == other)

    def doKwargsSubstitution(self, kwargs):
        oldStoichiometry = self.stoichiometry
        self.stoichiometry = {}
        for base in oldStoichiometry:
            self.stoichiometry[kwargs[base]] = oldStoichiometry[base]

        self.kineticLaw = ExprManip.sub_for_vars(self.kineticLaw, kwargs)

    def change_stoichiometry(self, species, stoich):
        """
        Change stoichiometry will update self.stoichiometry.
        It also updates self.reactant_stoichiometry and self.product_stoichiometry,
        if they are defined, to keep the two stoichiometry tracking systems in sync.
        """
        self.stoichiometry[species] = stoich
        if self.reactant_stoichiometry != None and self.product_stoichiometry != None:
            if species in self.reactant_stoichiometry.keys():
                del self.reactant_stoichiometry[species]
            self.product_stoichiometry[species] = [stoich]


class HomodimerizationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * reactant**2'
        self.stoichiometry = {'reactant': -2,
                              'dimer': +1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class HeterodimerizationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * A * B'
        self.stoichiometry = {'A': -1,
                              'B': -1,
                              'dimer': +1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class HomodimerDissociationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * dimer'
        self.stoichiometry = {'reactant': +2,
                              'dimer': -1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class HeterodimerDissociationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * dimer'
        self.stoichiometry = {'A': +1,
                              'B': +1,
                              'dimer': -1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class ExponentialDecayReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * species'
        self.stoichiometry = {'species': -1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class MichaelisMentenReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'k * E * S / (S + Km)'
        self.stoichiometry = {'E': 0,
                              'S': -1,
                              'P': 1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class ConstructionReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * template'
        self.stoichiometry = {'product': 1,
                              'template': 0} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class TransformationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * old'
        self.stoichiometry = {'old': -1,
                           'new': +1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class ProductionReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate'
        self.stoichiometry = {'product': 1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class PromoterReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'vmax * kP * P / (kP * P + ((kR1 * R1)**h) + ((kR2 * R2)**h) + 1.0)'
        self.stoichiometry = {'P': 0,
                              'R1': 0,
                              'R2': 0,
                              'mRNA': +1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class CoPromoterReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'vmax * (kP1 + kP2_1 * P) / (kP1 + kP2_2 * P + ((kR1 * R1 / (1 + kR2 * R2))**h) + 1.0)'
        self.stoichiometry = {'P': 0,
                              'R1': 0,
                              'R2': 0,
                              'mRNA': +1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class MichaelisMentenDegradationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'k * S / (S + Km)'
        self.stoichiometry = {'S': -1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class SelfCatalyticMichaelisMentenReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'k * S / (S + Km)'
        self.stoichiometry = {'S': -1,
                              'P': +1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class FirstOrderReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * reactant'
        self.stoichiometry = {'product': 1,
                              'reactant': -1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class TwoProductFirstOrderReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * reactant'
        self.stoichiometry = {'product_1': 1,
                              'product_2': 1,
                              'reactant': -1} 
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class SecondOrderReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'rate * reactant_1 * reactant_2'
        self.stoichiometry = {'product': 1,
                              'reactant_1': -1,
                              'reactant_2': -1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class HillDegradationReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'vmax * (S**h) / ((Km**h) + (S**h))'
        self.stoichiometry = {'S': -1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

class HillTransportReaction(Reaction):
    def __init__(self, id, **kwargs):
        self.kineticLaw = 'vmax * (S**h) / ((Km**h) + (S**h))'
        self.stoichiometry = {'S': -1,
                              'P': +1}
        self.doKwargsSubstitution(kwargs)
        Reaction.__init__(self, id, self.stoichiometry, self.kineticLaw)

