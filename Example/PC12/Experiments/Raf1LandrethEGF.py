import SloppyCell.Collections as Collections

expt = Collections.Experiment('Raf1LandrethEGF')

expt.longname = 'EGF Stimulation 100 ng/ml - Kao 2001'
expt.comments = """REF: S. Kao et. al., JBC, vol. 276, 18169 (2001)
CELLTYPE: PC12
MEAS: Raf1 activation in the presence of EGF at 100 ng/ml
UNITS: Relative kinase activity
NOTES: No error bars were quoted for this data set in the original paper
NOTES: 1.0 has been subtracted from this data so C(0) = 0 rather than 1.0"""

expt.SetData({'EGFstim100':{
                            'Raf1Active': {
                                           2.0:(2.75, 0.750),
                                           5.0:(1.6, 0.52),
                                           15.0:(0.25, 0.26),
                                           30.0:(0.0, 0.2),
                                           60.0:(0.0, 0.2),
                                           }
                            }
              }
             )
