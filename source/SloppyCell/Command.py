"""
Created 5/24/2017

Author @Keeyan <keeyan.ghoreshi@uconn.edu>

Controls command line options, input acceptance and output
"""

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
    root.update()   # First we update then we destroy
    root.destroy()  # It's what the documentation said to do don't ask questions
    parameters_to_fit, key_list = checkout(params, 'Select Parameters to Fit')  # list of user-chosen parameters
    for key in key_list:
        params.remove_by_key(key)
    parameters_to_skip, keys = checkout(params, 'Select Parameters to Skip')
    composition_dict['parameters_to_fit'] = parameters_to_fit
    composition_dict['parameters_to_skip'] = parameters_to_skip
    ScellParser.write_to_file(composition_dict)

# Main body start
parser = argparse.ArgumentParser(description='Accept SloppyCell inputs')
# Define parser arguments.
parser.add_argument('-i', type=str, dest='input_file', action='store', help='Specify input file location',
                    nargs='?', const='unspecified')
parser.add_argument('-o', type=str, dest='output_file', action='store', help='Specify output directory location')
parser.add_argument('--create', dest='file_creation', help='Create an .scell file.  Overrides other options.',
                    nargs='?', const=True, default=False)

input_v = parser.parse_args().input_file
output_v = parser.parse_args().output_file
file_create = parser.parse_args().file_creation
if file_create:
    create_scell()
