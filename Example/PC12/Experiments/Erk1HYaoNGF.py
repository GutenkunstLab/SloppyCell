import SloppyCell.Collections as Collections

expt = Collections.Experiment('Erk1HYaoNGF')

expt.longname = 'NGF Stimulation 50 ng/ml - Yao 1995'
expt.comments = """REF: H. Yao et. al., JBC 270(35), 20748
CELLTYPE: PC12
MEAS: Erk1 activation in the presence of NGF at 50 ng/ml
UNITS: Fold activation - 1 (so C(0) = 0)
NOTES: Error bars come from the original data"""

# XXX Kevin has two fixed scale factors for this expt, but no data for them?

expt.SetData({'NGFstim50': {
                            'ErkActive': {
                                          10.0:(6.8, 0.78),
                                          20.0:(7.8, 1.2),
                                          30.0:(7.1, 0.81),
                                          60.0:(6.6, 1.2)
                                          }
                            }
              }
             )
