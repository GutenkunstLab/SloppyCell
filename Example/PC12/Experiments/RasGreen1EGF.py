import SloppyCell.Collections as Collections

expt = Collections.Experiment('RasGreen1EGF')

expt.longname = 'EGF Stimulation 30 ng/ml - Qiu and Green 1994'
expt.comments = """REF: M-S. Qiu and S. H. Green, Neuron (1991) 7, 937-946
CELLTYPE: PC12
MEAS: GTP binding by Ras in response to EGF at 30 ng/ml
UNITS: Ras.GTP/(Ras.GTP + Ras.GDP) (percent total) X 10
NOTES: Error bars come from the original data set"""

expt.SetFixedScaleFactors({'RasActive': 1.0/(0.2*600000*0.1)})

expt.SetData({'EGFstim30': {
                            'RasActive': {
                                          3.0:(1.6, 0.4),
                                          10.0:(0.8, 	0.08),
                                          30.0:(0.45, 0.05)
                                          }
                            }
              }
             )
