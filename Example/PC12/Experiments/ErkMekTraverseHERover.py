import SloppyCell.Collections as Collections

expt = Collections.Experiment('ErkMekTraverseHERover')

expt.longname = '50-fold overexpression of EGFR, EGF 100 ng/ml - Traverse 1994'
expt.comments = """REF: S. Traverse, et. al., Curr. Biol. (1994) 4, 694
CELLTYPE: PC12
MEAS: Erk1/Mek activation in the presence of EGF at 100 ng/ml and
MEAS:    ~50-fold overexpression of wt HER (Human EGFR)
UNITS: units/mg (see paper), different def'n for Mek and Erk
NOTES: Error bars come directly from the data"""

expt.SetData({'EGFRx50_EGFstim100': {
                                     'ErkActive': {
                                                   2.0:(1.25, 0.25),
                                                   5.0:(3.5, 0.35),
                                                   15.0:(3.3, 0.35),
                                                   30.0:(2.5, 0.25),
                                                   60.0:(2.0, 0.25),
                                                   90.0:(1.4, 0.14)
                                                   },
                                     'MekActive': {
                                                   2.0:(0.9, 0.1),
                                                   5.0:(2.25, 0.25),
                                                   15.0:(0.8, 0.1),
                                                   30.0:(0.4, 0.1),
                                                   60.0:(0.35, 0.1),
                                                   90.0:(0.375, 0.05)
                                                   },
                                     }
              }
             )
