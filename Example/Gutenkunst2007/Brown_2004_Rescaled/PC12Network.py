import SloppyCell.ReactionNetworks
import SloppyCell.ReactionNetworks.Reactions as Reactions

baseNetwork = SloppyCell.ReactionNetworks.Network('base')

baseNetwork.addCompartment('cell')

baseNetwork.addSpecies('EGF', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('NGF', 'cell', 0.0, typicalValue = 1)
baseNetwork.addParameter('EGFR_IC', 1, isOptimizable=False)
baseNetwork.addSpecies('totalEGFReceptor', 'cell', 'EGFR_IC', isConstant=True)
baseNetwork.addSpecies('boundEGFReceptor', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('freeNGFReceptor', 'cell', 1)
baseNetwork.addSpecies('boundNGFReceptor', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('SosInactive', 'cell', 1)
baseNetwork.addSpecies('SosActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('P90RskInactive', 'cell', 1)
baseNetwork.addSpecies('P90RskActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('RasInactive', 'cell', 1)
baseNetwork.addSpecies('RasActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('RasGapActive', 'cell', 1, isConstant=True)
baseNetwork.addSpecies('Raf1Inactive', 'cell', 1)
baseNetwork.addSpecies('Raf1Active', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('BRafInactive', 'cell', 1)
baseNetwork.addSpecies('BRafActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('MekInactive', 'cell', 1)
baseNetwork.addSpecies('MekActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('ErkInactive', 'cell', 1)
baseNetwork.addSpecies('ErkActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('PI3KInactive', 'cell', 1)
baseNetwork.addSpecies('PI3KActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('AktInactive', 'cell', 1)
baseNetwork.addSpecies('AktActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('C3GInactive', 'cell', 1)
baseNetwork.addSpecies('C3GActive', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('Rap1Inactive', 'cell', 1)
baseNetwork.addSpecies('Rap1Active', 'cell', 0.0, typicalValue = 1)
baseNetwork.addSpecies('RapGapActive', 'cell', 1, isConstant = True)
baseNetwork.addSpecies('PP2AActive', 'cell', 1, isConstant = True)
baseNetwork.addSpecies('Raf1PPtase', 'cell', 1, isConstant=True)

baseNetwork.addParameter('kdEGF', 0, isOptimizable=False)
baseNetwork.addParameter('krbNGF', 0.07673)
baseNetwork.addParameter('kruNGF', 0, isOptimizable=False)

baseNetwork.addParameter('kEGF', 10.67)
baseNetwork.addParameter('kNGF', 14.43)

baseNetwork.addParameter('kdSos', 21.35)
baseNetwork.addParameter('kSos', 0.7701)

baseNetwork.addParameter('kRasGap', 2.310)
baseNetwork.addParameter('kRasToRaf1', 2.264)

baseNetwork.addParameter('kpRaf1', 1.633)
baseNetwork.addParameter('kpBRaf', 0.1015)

baseNetwork.addParameter('kdMek', 0.6177)

baseNetwork.addParameter('kpMekCytoplasmic', 21.35)
baseNetwork.addParameter('kdErk', 0.5011)
baseNetwork.addParameter('kpP90Rsk', 0.2716)

baseNetwork.addParameter('kPI3K', 21.35)
baseNetwork.addParameter('kPI3KRas', 1.660)
baseNetwork.addParameter('kAkt', 0.07673)
baseNetwork.addParameter('kdRaf1ByAkt', 13.32)

baseNetwork.addParameter('kC3GNGF', 3.847)
baseNetwork.addParameter('kC3G', 0.1976)
baseNetwork.addParameter('kRapGap', 0.8229)
baseNetwork.addParameter('kRap1ToBRaf', 21.35)

baseNetwork.addParameter('kdRaf1', 1.557)

baseNetwork.addParameter('kdBRaf', 0.4363)

# Put EGF receptor in equilibrium
baseNetwork.addAssignmentRule('boundEGFReceptor', '0.5*(EGF+totalEGFReceptor+kdEGF - sqrt((EGF + totalEGFReceptor + kdEGF)**2 - 4*EGF*totalEGFReceptor))')

# At least for the best-fit parameters, it appears necessary to have NGF not
#  immediately equilibrate.
#baseNetwork.addParameter('kdNGF', isOptimizable=False)
#baseNetwork.addAssignmentRule('kdNGF', 'kruNGF/krbNGF')
#baseNetwork.addAssignmentRule('boundNGFReceptor', '0.5*(NGF+freeNGFReceptor+kdNGF - sqrt((NGF + freeNGFReceptor + kdNGF)**2 - 4*NGF*freeNGFReceptor))')
baseNetwork.addReaction(Reactions.HeterodimerizationReaction,
                        'NGFBindingReaction',
                        A='NGF', B='freeNGFReceptor', dimer='boundNGFReceptor',
                        rate='krbNGF')

baseNetwork.addReaction(Reactions.HeterodimerDissociationReaction,
                        'NGFUnbindingReaction',
                        dimer='boundNGFReceptor', A='freeNGFReceptor', B='NGF',
                        rate='kruNGF')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'SosActivationByEGFReaction',
                        E='boundEGFReceptor', S='SosInactive', P='SosActive',
                        k='kEGF', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'SosActivationByNGFReaction',
                        E='boundNGFReceptor', S='SosInactive', P='SosActive',
                        k='kNGF', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'SosDeactivationReaction',
                        E='P90RskActive', S='SosActive', P='SosInactive',
                        k='kdSos', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'RasActivationReaction',
                        E='SosActive', S='RasInactive', P='RasActive',
                        k='kSos', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'RasDeactivationReaction',
                        E='RasGapActive', S='RasActive', P='RasInactive',
                        k='kRasGap', Km=1)

###
### RNG7: Modified activation of Raf1 and PI3K by Ras to account for fact
###       that there are *two* substrates for one enzyme.
###
ratelaw = '%(k1)s*%(E)s * %(S1)s * (1-%(S2)s/(%(kd2)s + %(S2)s))/(%(kd1)s + %(S1)s*(1 + %(S2)s/(%(kd2)s + %(S2)s)))'
baseNetwork.addReaction('Raf1ByRasActivationReaction',
                        {'RasActive':0, 'Raf1Inactive':-1, 'Raf1Active':1},
                        ratelaw % {'E':'RasActive',
                                   'k1': 'kRasToRaf1',
                                   'S1': 'Raf1Inactive',
                                   'kd1': 1,
                                   'S2': 'PI3KInactive',
                                   'kd2': 1})
baseNetwork.addReaction('PI3KByRasActivationReaction',
                        {'RasActive':0, 'PI3KInactive':-1, 'PI3KActive':1},
                        ratelaw % {'E':'RasActive',
                                   'k1': 'kPI3KRas',
                                   'S1': 'PI3KInactive',
                                   'kd1': 1,
                                   'S2': 'Raf1Inactive',
                                   'kd2': 1,})


baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekbyRaf1ActivationReaction',
                        E='Raf1Active', S='MekInactive', P='MekActive',
                        k='kpRaf1', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekbyBRafActivationReaction',
                        E='BRafActive', S='MekInactive', P='MekActive',
                        k='kpBRaf', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'ErkActivationReaction',
                        E='MekActive', S='ErkInactive', P='ErkActive',
                        k='kpMekCytoplasmic', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekDeactivationReaction',
                        E='PP2AActive', S='MekActive', P='MekInactive',
                        k='kdMek', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'ErkDeactivationReaction',
                        E='PP2AActive', S='ErkActive', P='ErkInactive',
                        k='kdErk', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Raf1byPPtaseDeactivationReaction',
                        E='Raf1PPtase', S='Raf1Active', P='Raf1Inactive',
                        k='kdRaf1', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'P90RskActivationReaction',
                        E='ErkActive', S='P90RskInactive', P='P90RskActive',
                        k='kpP90Rsk', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'PI3KbyEGFRActivationReaction',
                        E='boundEGFReceptor', S='PI3KInactive', P='PI3KActive',
                        k='kPI3K', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'AktActivationReaction',
                        E='PI3KActive', S='AktInactive', P='AktActive',
                        k='kAkt', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Raf1ByAktDeactivationReaction',
                        E='AktActive', S='Raf1Active', P='Raf1Inactive',
                        k='kdRaf1ByAkt', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'C3GActivationReaction',
                        E='boundNGFReceptor', S='C3GInactive', P='C3GActive',
                        k='kC3GNGF', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Rap1ActivationReaction',
                        E='C3GActive', S='Rap1Inactive', P='Rap1Active',
                        k='kC3G', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Rap1DeactivationReaction',
                        E='RapGapActive', S='Rap1Active', P='Rap1Inactive',
                        k='kRapGap', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'BRafByRap1ActivationReaction',
                        E='Rap1Active', S='BRafInactive', P='BRafActive',
                        k='kRap1ToBRaf', Km=1)

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'BRafbyPPtaseDeactivationReaction',
                        E='Raf1PPtase', S='BRafActive', P='BRafInactive',
                        k='kdBRaf', Km=1)

baseNetwork.compile()
