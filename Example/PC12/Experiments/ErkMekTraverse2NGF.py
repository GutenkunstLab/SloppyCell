import SloppyCell.Collections as Collections

expt = Collections.Experiment('ErkMekTraverse2NGF')

expt.longname = 'NGF Stimulation 50 ng/ml - Traverse 1994'
expt.comments = """REF: S. Traverse, et. al., Curr. Biol. (1994) 4, 694
CELLTYPE: PC12
MEAS: Erk1/Mek activation in the presence of NGF at 50 ng/ml
UNITS: units/mg (see paper), different def'n for Erk and Mek
NOTES: Error bars are those in the original data"""

expt.SetData({'NGFstim50': {
                            'MekActive': {
                                          2.0:(0.375, 0.05),
                                          5.0:(1.8, 0.3),
                                          15.0:(1.13, 0.113),
                                          30.0:(1.4, 0.14),
                                          60.0:(1.3, 0.13),
                                          90.0:(0.8, 0.13)
                                          },
                            'ErkActive': {
                                          2.0:(0.4, 0.1),
                                          5.0:(3.7, 0.37),
                                          15.0:(3.5, 0.4),
                                          30.0:(4.1, 0.41),
                                          60.0:(3.75, 0.5),
                                          90.0:(3.1, 0.31)
                                          }
                            }
              }
             )
