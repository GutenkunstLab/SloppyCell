import SloppyCell.Collections as Collections

expt = Collections.Experiment('ErkMekTraverse2EGF')

expt.longname = 'EGF Stimulation 100 ng/ml - Traverse 1994'
expt.comments = """REF: S. Traverse, et. al., Curr. Biol. (1994) 4, 694
CELLTYPE: PC12
MEAS: Erk1/Mek activation in the presence of EGF at 100 ng/ml
UNITS: units/mg (see paper), different def'n for Mek and Erk
NOTES: Error bars are those in the original data, not made up
NOTES:"""

expt.SetData({'EGFstim100':{
                            'ErkActive': {
                                          2.0:(2.85, 0.29),
                                          5.0:(4.9, 0.49),
                                          15.0:(2.35, 0.23),
                                          30.0:(1.9, 0.35),
                                          60.0:(0.9, 0.10),
                                          90.0:(0.5, 0.10)
                                          },
                            'MekActive': {
                                          2.0:(1.25, 0.13),
                                          5.0:(1.12, 0.13),
                                          15.0:(0.2, 0.05),
                                          30.0:(0.18, 0.05),
                                          60.0:(0.125, 0.05),
                                          90.0:(0.1, 0.05)
                                          },
                            }
              }
             )
