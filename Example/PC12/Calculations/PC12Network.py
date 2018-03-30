import SloppyCell.ReactionNetworks
import SloppyCell.ReactionNetworks.Reactions as Reactions

CFactor = 600000.0
SecToMin = 60.0

baseNetwork = SloppyCell.ReactionNetworks.Network('base')

baseNetwork.addCompartment('cell')

baseNetwork.addSpecies('EGF', 'cell', 0.0, typicalValue = 1e8)
baseNetwork.addSpecies('NGF', 'cell', 0.0, typicalValue = 1e6)
baseNetwork.addSpecies('freeEGFReceptor', 'cell', 80000.0)
baseNetwork.addSpecies('boundEGFReceptor', 'cell', 0.0, typicalValue = 1e4)
baseNetwork.addSpecies('freeNGFReceptor', 'cell', 10000.0)
baseNetwork.addSpecies('boundNGFReceptor', 'cell', 0.0, typicalValue = 1e4)
baseNetwork.addSpecies('SosInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('SosActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('P90RskInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('P90RskActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('RasInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('RasActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('RasGapActive', 'cell', 0.2*CFactor)#, isConstant=True)
baseNetwork.addSpecies('Raf1Inactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('Raf1Active', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('BRafInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('BRafActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('MekInactive', 'cell', CFactor)
baseNetwork.addSpecies('MekActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('ErkInactive', 'cell', CFactor)
baseNetwork.addSpecies('ErkActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('PI3KInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('PI3KActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('AktInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('AktActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('C3GInactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('C3GActive', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('Rap1Inactive', 'cell', 0.2*CFactor)
baseNetwork.addSpecies('Rap1Active', 'cell', 0.0, typicalValue = CFactor)
baseNetwork.addSpecies('RapGapActive', 'cell', 0.2*CFactor)#, isConstant = True)
baseNetwork.addSpecies('PP2AActive', 'cell', 0.2*CFactor)#, isConstant = True)
baseNetwork.addSpecies('Raf1PPtase', 'cell', 0.2*CFactor)#, isConstant=True)

baseNetwork.addParameter('krbEGF', SecToMin*0.00001)
baseNetwork.addParameter('kruEGF', SecToMin*0.01)
baseNetwork.addParameter('krbNGF', SecToMin*0.00005)
baseNetwork.addParameter('kruNGF', SecToMin*0.02)

baseNetwork.addParameter('kEGF', SecToMin*0.2)
baseNetwork.addParameter('KmEGF', CFactor*0.8333)
baseNetwork.addParameter('kNGF', SecToMin*0.4)
baseNetwork.addParameter('KmNGF', CFactor*0.8333)

baseNetwork.addParameter('kdSos', 0.01*SecToMin)
baseNetwork.addParameter('KmdSos', 100000.0)
baseNetwork.addParameter('kSos', 0.075*SecToMin)
baseNetwork.addParameter('KmSos', 0.25*CFactor)

baseNetwork.addParameter('kRasGap', 0.2*SecToMin)
baseNetwork.addParameter('KmRasGap', CFactor*1.0104)
baseNetwork.addParameter('kRasToRaf1', 0.08*SecToMin)
baseNetwork.addParameter('KmRasToRaf1', 150000.0)

baseNetwork.addParameter('kpRaf1', 0.05*SecToMin)
baseNetwork.addParameter('KmpRaf1', 0.6*CFactor)
baseNetwork.addParameter('kpBRaf', 0.025*SecToMin)
baseNetwork.addParameter('KmpBRaf', 100000.0)

baseNetwork.addParameter('kdMek', 0.03*SecToMin)
baseNetwork.addParameter('KmdMek', 125000.0)

baseNetwork.addParameter('kpMekCytoplasmic', 0.01*SecToMin)
baseNetwork.addParameter('KmpMekCytoplasmic', 125000.0)
baseNetwork.addParameter('kdErk', 0.01*SecToMin)
baseNetwork.addParameter('KmdErk', 300000.0)
baseNetwork.addParameter('kpP90Rsk', 0.01*SecToMin)
baseNetwork.addParameter('KmpP90Rsk', 300000.0)

baseNetwork.addParameter('kPI3K', 0.05*SecToMin)
baseNetwork.addParameter('KmPI3K', 100000.0)
baseNetwork.addParameter('kPI3KRas', 0.025*SecToMin)
baseNetwork.addParameter('KmPI3KRas', 300000.0)
baseNetwork.addParameter('kAkt', 0.01*SecToMin)
baseNetwork.addParameter('KmAkt', 300000.0)
baseNetwork.addParameter('kdRaf1ByAkt', 0.1*SecToMin)
baseNetwork.addParameter('KmRaf1ByAkt', 300000.0)

baseNetwork.addParameter('kC3GNGF', 0.05*SecToMin)
baseNetwork.addParameter('KmC3GNGF', 200000.0)
baseNetwork.addParameter('kC3G', 0.01*SecToMin)
baseNetwork.addParameter('KmC3G', 300000.0)
baseNetwork.addParameter('kRapGap', 0.002*SecToMin)
baseNetwork.addParameter('KmRapGap', 300000.0)
baseNetwork.addParameter('kRap1ToBRaf', 0.05*SecToMin)
baseNetwork.addParameter('KmRap1ToBRaf', 300000.0)

baseNetwork.addParameter('kdRaf1', 0.1*SecToMin)
baseNetwork.addParameter('KmdRaf1', 300000.0)

baseNetwork.addParameter('kdBRaf', 0.05*SecToMin)
baseNetwork.addParameter('KmdBRaf', 300000.0)

baseNetwork.addReaction(Reactions.HeterodimerizationReaction,
                        'EGFBindingReaction',
                        A='EGF', B='freeEGFReceptor', dimer='boundEGFReceptor',
                        rate='krbEGF')

baseNetwork.addReaction(Reactions.HeterodimerDissociationReaction,
                        'EGFUnbindingReaction',
                        dimer='boundEGFReceptor', A='freeEGFReceptor', B='EGF',
                        rate='kruEGF')

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
                        k='kEGF', Km='KmEGF')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'SosActivationByNGFReaction',
                        E='boundNGFReceptor', S='SosInactive', P='SosActive',
                        k='kNGF', Km='KmNGF')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'SosDeactivationReaction',
                        E='P90RskActive', S='SosActive', P='SosInactive',
                        k='kdSos', Km='KmdSos')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'RasActivationReaction',
                        E='SosActive', S='RasInactive', P='RasActive',
                        k='kSos', Km='KmSos')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'RasDeactivationReaction',
                        E='RasGapActive', S='RasActive', P='RasInactive',
                        k='kRasGap', Km='KmRasGap')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Raf1ByRasActivationReaction',
                        E='RasActive', S='Raf1Inactive', P='Raf1Active',
                        k='kRasToRaf1', Km='KmRasToRaf1')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekbyRaf1ActivationReaction',
                        E='Raf1Active', S='MekInactive', P='MekActive',
                        k='kpRaf1', Km='KmpRaf1')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekbyBRafActivationReaction',
                        E='BRafActive', S='MekInactive', P='MekActive',
                        k='kpBRaf', Km='KmpBRaf')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'ErkActivationReaction',
                        E='MekActive', S='ErkInactive', P='ErkActive',
                        k='kpMekCytoplasmic', Km='KmpMekCytoplasmic')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'MekDeactivationReaction',
                        E='PP2AActive', S='MekActive', P='MekInactive',
                        k='kdMek', Km='KmdMek')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'ErkDeactivationReaction',
                        E='PP2AActive', S='ErkActive', P='ErkInactive',
                        k='kdErk', Km='KmdErk')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Raf1byPPtaseDeactivationReaction',
                        E='Raf1PPtase', S='Raf1Active', P='Raf1Inactive',
                        k='kdRaf1', Km='KmdRaf1')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'P90RskActivationReaction',
                        E='ErkActive', S='P90RskInactive', P='P90RskActive',
                        k='kpP90Rsk', Km='KmpP90Rsk')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'PI3KbyEGFRActivationReaction',
                        E='boundEGFReceptor', S='PI3KInactive', P='PI3KActive',
                        k='kPI3K', Km='KmPI3K')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'PI3KbyRasActivationReaction',
                        E='RasActive', S='PI3KInactive', P='PI3KActive',
                        k='kPI3KRas', Km='KmPI3KRas')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'AktActivationReaction',
                        E='PI3KActive', S='AktInactive', P='AktActive',
                        k='kAkt', Km='KmAkt')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Raf1ByAktDeactivationReaction',
                        E='AktActive', S='Raf1Active', P='Raf1Inactive',
                        k='kdRaf1ByAkt', Km='KmRaf1ByAkt')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'C3GActivationReaction',
                        E='boundNGFReceptor', S='C3GInactive', P='C3GActive',
                        k='kC3GNGF', Km='KmC3GNGF')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Rap1ActivationReaction',
                        E='C3GActive', S='Rap1Inactive', P='Rap1Active',
                        k='kC3G', Km='KmC3G')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'Rap1DeactivationReaction',
                        E='RapGapActive', S='Rap1Active', P='Rap1Inactive',
                        k='kRapGap', Km='KmRapGap')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'BRafByRap1ActivationReaction',
                        E='Rap1Active', S='BRafInactive', P='BRafActive',
                        k='kRap1ToBRaf', Km='KmRap1ToBRaf')

baseNetwork.addReaction(Reactions.MichaelisMentenReaction,
                        'BRafbyPPtaseDeactivationReaction',
                        E='Raf1PPtase', S='BRafActive', P='BRafInactive',
                        k='kdBRaf', Km='KmdBRaf')

baseNetwork.compile()
