import SloppyCell.Collections as Collections

expt = Collections.Experiment('Raf1BrafLandrethNGF')

expt.longname = 'NGF Stimulation 100 ng/ml - Kao 2001'
expt.comments = """REF: S. Kao et. al., JBC, vol. 276, 18169 (2001)
CELLTYPE: PC12
MEAS: Raf1 activation in the presence of NGF at 100 ng/ml
MEAS: BRaf activation in the presence of NGF at 100 ng/ml
UNITS: Relative kinase activity
NOTES: Error bars are read from the figures as accurately as I am able.
NOTES: 1.0 has been subtracted from this data so C(0) = 0 rather than 1.0"""

expt.SetData({'NGFstim100': {
                             'Raf1Active': {
                                            2.0:(1.8, 0.28),
                                            5.0:(2.6, 0.36),
                                            15.0:(0.6, 0.16),
                                            30.0:(0.25, 0.25),
                                            60.0:(0.1, 0.20)
                                            },
                             'BRafActive': {
                                            2.0:(1.25, 0.22),
                                            5.0:(1.9, 0.29),
                                            15.0:(1.4, 0.24),
                                            30.0:(1.2, 0.25),
                                            60.0:(1.1, 0.21)
                                            }
                             }
              }
             )
