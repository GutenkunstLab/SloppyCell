import SloppyCell.Collections as Collections

expt = Collections.Experiment('RasGreen1NGF')

expt.longname = 'NGF Stimulation 50 ng/ml - Qiu and Green 1994'
expt.comments = """REF: M-S. Qiu and S. H. Green, Neuron (1991) 7, 937-946
CELLTYPE: PC12
MEAS: GTP binding by Ras in response to NGF at 50 ng/ml
UNITS: Ras.GTP/(Ras.GTP + Ras.GDP) X 10 (percent total / 10)
NOTES: Multiplication by 10 is just to put data on a better scale for
NOTES:    plotting
NOTES: Error bars come from the original data set"""

expt.SetFixedScaleFactors({'RasActive':1.0/(0.2*600000*0.1)})

expt.SetData({'NGFstim50': {
                            'RasActive': {
                                          3.0:(1.7, 0.5),
                                          10.0:(1.1, 0.2),
                                          30.0:(0.6, 0.1)
                                          }
                            }
              }
             )
