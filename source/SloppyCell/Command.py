"""
Created 5/24/2017

Author @Keeyan <keeyan.ghoreshi@uconn.edu>

Controls command line options, input acceptance and output
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # Adds root directory to path
import argparse
from Tkinter import *
from tkFileDialog import *
import ScellParser
import ReactionNetworks as ReactionNetworks
import logging
logger = logging.getLogger('Command')


class Checkbar(Frame):
    # Class for extensible bar of checkboxes (Tkinter)
    def __init__(self, parent=None, picks=[], side=TOP, anchor=W):
        Frame.__init__(self, parent)
        self.vars = []
        for pick in picks:
            var = IntVar()
            chk = Checkbutton(self, text=pick, variable=var)
            chk.pack(side=side, anchor=anchor, expand=YES)
            self.vars.append(var)

    def state(self):
        return map((lambda var: var.get()), self.vars)


def checkout(params, title):
    """
    Accepts a list of parameters and a window title.
    Creates a window with the chosen title and a column of 
    checkboxes next to each parameter in the list.
    
    Returns two lists, one of the keys (or parameter names) and another
    of all the values of each parameter.
    """
    root = Tk()
    root.attributes("-topmost", True)  # bring the window to front
    root.lift()
    w = Label(root, text=title, bd=3)
    w.pack(side=TOP)
    lng = Checkbar(root, params.keys())

    lng.pack(side=TOP,  fill=Y)
    lng.config(relief=GROOVE, bd=2)

    Button(root, text='Select Parameters', command=root.quit, bd=2,
           default=ACTIVE).pack(side=BOTTOM, fill=BOTH, expand=YES)

    # Button(root, text='Peek', command=allstates).pack(side=RIGHT)  # Useful for debugging
    root.mainloop()  # Program halts until window is closed

    input_list = lng.state()
    output_list = []
    key_list = []
    for x in range(len(input_list)):
        if list(input_list)[x] == 1:
            output_list.append(params.items()[x])
            key_list.append(params.keys()[x])

    root.update()  # First we update then we destroy
    root.destroy()
    return output_list, key_list

def create_scell():
    """
    Returns a dictionary that can be used
    to build a text file of the format ->
    <Category>[value]
    """
    # TODO: update creation function to include new features.  Currently: initial condition setting, prior bounds
    composition_dict = {}
    root = Tk()  # Create window for dialog boxes
    root.withdraw()  # hide root window
    root.update()
    sbml_filename = askopenfilename(title='Select an SBML file')  # show dialog box for directory selection
    composition_dict['SBML_Reference'] = sbml_filename
    root.update()
    data_filename = askopenfilename(title='Select a data file')
    composition_dict['Data_Reference'] = data_filename
    net = ReactionNetworks.IO.from_SBML_file(sbml_filename)
    params = net.GetParameters()
    print net.parameters
    root.update()   # First we update then we destroy
    root.destroy()  # It's what the documentation said to do don't ask questions
    parameters_to_fit, key_list = checkout(params, 'Select Parameters to Fit')  # list of user-chosen parameters
    for key in key_list:
        params.remove_by_key(key)
    parameters_to_skip, keys = checkout(params, 'Select Parameters to Skip')
    composition_dict['parameters_to_fit'] = parameters_to_fit
    composition_dict['parameters_to_skip'] = parameters_to_skip
    return composition_dict

# Main body start
# should probably be inside a main function
parser = argparse.ArgumentParser(description='Accept SloppyCell inputs')
# Define parser arguments.
# Quick explanation:
# first argument is the flag
# >type< is type of input accepted
# >dest< is the name of the variable used to access arguments put after flag
# >action< is what is to be done with the arguments put after flag (usually just store)
# >help< is the text that will be displayed when the --help command is typed
# >nargs< is how many arguments are expected, and whether the command is required ('?' means that flag is not required)
# >const< is the value the flag returns if it is specified but nothing is put after it (no arguments given)
# >default< is the value the flag returns if it is not specified at all.
parser.add_argument('-i', type=str, dest='input_file', action='store', help='Specify input file location',
                    nargs='?', default='unspecified')
parser.add_argument('-o', type=str, dest='output_file', action='store', help='Specify output directory location')
parser.add_argument('--create', dest='file_creation', help='Create an .scell file.  Overrides other options.',
                    nargs='?', const=True, default=False)

input_v = parser.parse_args().input_file
output_v = parser.parse_args().output_file
file_create = parser.parse_args().file_creation
if file_create:
    ScellParser.write_to_file(create_scell())
if input_v is not 'unspecified':
    ScellParser.read_from_file(input_v)
