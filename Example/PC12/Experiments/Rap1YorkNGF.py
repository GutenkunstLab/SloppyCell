import SloppyCell.Collections as Collections

expt = Collections.Experiment('Rap1YorkNGF')

expt.longname = 'NGF Stimulation 50 ng/ml - York 1998'
expt.comments = """REF: R. D. York et. al., Nature vol. 392, 622 (1998)
CELLTYPE: PC12
MEAS: GTP binding by Rap1 in response to NGF at 50 ng/ml
UNITS: Rap1.GTP/(Rap1.GTP + Rap1.GDP) X 10 (percent total / 10)
NOTES: Division by 10 is just to put data on a better scale for
NOTES:    plotting
NOTES: Rap1(0) = 0.48 subtracted from this data set so that Rap1(0) = 0
NOTES: No error bars were provided with the original data"""

expt.SetFixedScaleFactors({'Rap1Active': 1.0/(0.2*600000*0.1)})

expt.SetData({'NGFstim50': {
                            'Rap1Active': {
                                           #   0.0:(Rap1Active    0.48   0.048
                                           10.0:(1.09, 0.314),
                                           30.0:(1.11, 0.318),
                                           40.0:(1.59, 0.414),
                                           }
                            }
              }
             )
