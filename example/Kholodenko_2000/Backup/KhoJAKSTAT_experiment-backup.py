from __future__ import print_function
import numpy as np
from numpy import genfromtxt
import csv
import sys

f = open('data.csv', 'r')
first_line = f.readline()
cells = first_line.split( "," )
num_var=len(cells)
var_count=0
print ("Number of vartiables", num_var, "\n")
f.close()

print ("Opening the file...\n")
target = open("C:\\Users\\vcguest\\Desktop\\Sloppy\\Sloppy Cell Download\\Example\\Kholodenko_2000\\output.txt", 'w')
print ("from SloppyCell.ReactionNetworks import *\n\nexpt = Experiment('expt1')\n", file=target)
for i in range (0, num_var-2):
    f = open('data.csv', 'r')
    var_count=var_count+1
    for line in f:
        cells = line.split( "," )
        if 'Time' in cells:
            print("data = {'net1':{'",cells[var_count], "': {    ", file=target)
        else:
            print("                      ", cells[0] ,": (",cells[var_count],", 0.1),", file=target)
        
    f.close()
    
print("                          }\n                }\n        }\nexpt.set_data(data)", file=target)
target.close()
        


    

    



