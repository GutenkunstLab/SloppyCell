"""
Created 5/24/2017

Author @Keeyan <keeyan.ghoreshi@uconn.edu>

Controls command line options, input acceptance and output
"""

import os
import sys


sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # Adds root directory to path
import argparse
import ScellParser
import logging
logger = logging.getLogger('Command')

# Type "python command.py -i <file_to_input>" into command line to run SloppyCell



# Main body start
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
print input_v
if input_v is not 'unspecified':
    ScellParser.read_from_file(input_v)
