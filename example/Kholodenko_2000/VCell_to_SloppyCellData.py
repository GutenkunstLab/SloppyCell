from __future__ import print_function
import numpy as np
import random
from numpy import genfromtxt
import csv

f = open('data.csv', 'r')
first_line = f.readline()
cells = first_line.split( "," )
num_var=len(cells);
num_times=0
for line in f:
        cells = line.split( "," )
        num_var1=len(cells);
        if num_var1 != num_var :
            print("Different number of columns in lines: exiting")
            exit()
        else:
            num_times+=1
f.close()
print("Number of time points is ",num_times )

var_count=0;
print ("Number of variables", num_var-1, "\n")
f.close()

print ("Opening the file...\n")
target = open("data_to_fit.py", 'w')

print ("from SloppyCell.ReactionNetworks import *\n\nexpt = Experiment('expt1')", file=target, sep='')
first_line_string = "data = {'net1':{'"
for i in range (0, num_var-1):
    f = open('data.csv', 'r')
    var_count=var_count+1
    line_count=0
    string = ""
    for line in f:
        cells = line.split( "," )
        if line_count > 0:
            variable1 = lambda: random.randint(-int(10e5*float(cells[var_count])), int(10e5*float(cells[var_count])))
            random_value = float(cells[var_count]) + variable1()/10e6 
        if line_count==0:
            string = cells[var_count].strip()
        elif line_count==1 and var_count == 1:
            print(first_line_string, string, "':{",cells[0] ,": (" , random_value, ", 0.1),", file=target, sep='')          
        elif  line_count==1:
            print("                '",string,"':{",cells[0] ,": (" , random_value,  ", 0.1),", file=target, sep='')
        elif  line_count==num_times:
            print("                      ",cells[0] ,": (" , random_value,", 0.1),\n                     },", file=target, sep='')
        elif line_count > 1:
            print("                      ",cells[0] ,": (" , random_value,", 0.1),", file=target, sep='')
        line_count+=1
    f.close()
print("                }\n        }\nexpt.set_data(data)", file=target)
target.close()


