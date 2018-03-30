import SloppyCell.Collections as Collections

expt = Collections.Experiment('Erk1HYaoEGF')

expt.longname = 'EGF Stimulation 100 ng/ml - Yao 1995'
expt.comments = """REF: H. Yao et. al., JBC 270(35), 20748
CELLTYPE: PC12
MEAS: Erk1 activation in the presence of EGF at 100 ng/ml
UNITS: Fold activation - 1 (so C(0) = 0), measured with a phosphorimager
NOTES: Error bars come from the original data"""

expt.SetData({'EGFstim100': {
                             "ErkActive": {
                                           10.0:(9.0, 1.0),
                                           20.0:(6.4, 0.74),
                                           30.0:(2.0, 0.3),
                                           60.0:(1.5, 0.3)
                                           }
                             }
              }
             )
