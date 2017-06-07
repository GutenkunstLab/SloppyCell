"""
Created 5/24/2017

Author @Keeyan

Parses .scell filetypes
"""
import os
import TestConstruct
import TestConstruct_XML
import logging
import re
import csv

from SloppyCell.ReactionNetworks import *
from xml.etree import ElementTree as ET


logger = logging.getLogger('ScellParser')
logging.basicConfig()
# TODO: Defining global variables is bad practice so change this later.  Should really just make it a class. @Keeyan
tag_string = '<>'  # change things here to reformat Scell files
data_tag_string = '[]'
accepted_data_file_type = ('.csv', 'tsv', '.xlsx', '.txt')
accepted_input_file_type = ('.xml', '.txt')


def tag(category, data):
    """
    Surrounds strings with the appropriate tags.
    Organizes data in a file to be read later.
    Returns two strings, a category and a data string.  Categories are contained 
    using the globally defined 'tag_string', and data are contained using the 
    'data_tag_string'.
    """
    return_string_d = ''  # must be mutable

    return_string_c = tag_string[0] + category + tag_string[1]
    try:
        return_string_d = data_tag_string[0] + data + data_tag_string[1]
    except TypeError:  # Data is sometimes a list of tuples with 2 values, handled here
        for pair in data:
            return_string_d += pair[0] + ',' + str(pair[1]) + ','
        excess = len(return_string_d)-2  # Index used to cut extra comma and period off end of each data string
        return_string_d = data_tag_string[0] + return_string_d[:excess] + data_tag_string[1]

    return return_string_c, return_string_d


def retuple(prelim_dict):
    """
    It's called retuple but it actually makes a dictionary (it's called retuple because parameters are originally
    held in a parameter:value tuple, but it makes more sense to put them in a dictionary within a dictionary)
    Finds 'parameter' lists and makes them their own dictionary
    """

    # Python is pass-by-object-reference, so modifying the object will work for references
    # elsewhere in the program, but reassignment will not register.
    # This could be a problem for maintenance/updates later, even if it works right now
    matching = [s for s in prelim_dict.keys() if "parameter" in s]
    beta = ''
    for key in matching:
        params = prelim_dict[key]
        param_list = re.split(',', params)
        placeholder = []
        for param in param_list:
            # This uses an exception to generate a parameter:value dictionary.  Probably not the best practice...
            # But it works.
            try:
                alpha = float(param)
                placeholder.append((beta, alpha))
            except ValueError:
                beta = param
        prelim_dict[key] = KeyedList(placeholder)


def untag(file_string):
    """
    Turns scell file contents into a dictionary. 
    
    Returns the dictionary it makes.
    
    
    """

    # TODO: get rid of regex functions
    return_dictionary = {}
    all_items = re.findall('<.*>\[.*\]', file_string)  # Makes a list of every item in the file
    # .* denotes any character, and any amount.  Essentially the regex expression checks for characters enclosed by
    # '<>' followed by '[]' (the brackets need to be escaped)
    for item in all_items:
        pair = re.split('>\[', item)
        return_dictionary[pair[0][1:]] = pair[1][:len(pair[1])-1]  # Chops off leading '<' and trailing ']'
    retuple(return_dictionary)
    return return_dictionary


def experiment_constructor(data_file):
    """
    Uses provided data file path to find appropriate csv file
    Extracts data and formats it to create an 'Experiment' object
    Returns the experiment object to be utilized in the main function
    """
    expt = Experiment(os.path.splitext(os.path.basename(data_file))[0])  # Experiment named after data file
    # expt = Experiment('expt1')
    if data_file.lower().endswith(accepted_data_file_type):
        with open(data_file) as csv_file:
            reader_v = csv.DictReader(csv_file)
            firstline = reader_v.next()
            keys = firstline.keys()
            # Read each line, categorize data by header into a dictionary
            model_dict = {}
            for key in keys:
                model_dict[key] = []

            # Deals with first line case
            for key in keys:
                model_dict[key].append(firstline[key])
            # Organizes lines into dictionaries, with keys referencing lists of all values associated with that key
            for row in reader_v:
                for key in keys:
                    model_dict[key].append(row[key])
            time_key = keys[len(keys)-1]  # Assumes that the TIME column is first, and a column
            keys.remove(keys[len(keys)-1])
            final_dict = {}
            time_array = []

            for key in keys:
                try:
                    data = model_dict[key]
                    temp_dict = {}
                    index = 0
                    for data_point in data:
                        # TODO: Hopefully find a way to avoid hard-coding in sigma-value (Technically done already)
                        sigma_list = data_point.split(',')
                        temp_dict[float(model_dict[time_key][index])] = (float(sigma_list[0]), float(sigma_list[1]))
                        index += 1
                    final_dict[key] = temp_dict
                    time_array = temp_dict.keys()

                except Exception as e:
                    print e
                    print 'This exception was brought to you by bad csv parsing'
        return_dict = dict()
        # TODO: Switch from hard-coded model inclusion to SBML file reading
        return_dict['net1'] = final_dict
        expt.set_data(return_dict)
        return expt, time_array
    else:
        logger.warn('Selected data file type not supported.')


def write_to_file(package):
    """
    Accepts a dictionary as a parameter.
    Writes it to a file with a 
    <Key>[Value] format.
    """
    target = open('testFile', 'w')
    for key in package.keys():
        a, b = tag(key, package[key])  # Adds the <> and [] tags around the key and value pair
        target.write(a + b + '\n')
    target.close()


def read_from_file(file_name):
    # try:
        if file_name.lower().endswith('.xml'):
            xml_file = ET.parse(file_name)
            root = xml_file.getroot()
            try:
                data_reference = root.find("References").find("Data").attrib["path"]
                experiment, time_array = experiment_constructor(data_reference)
            except Exception as e:
                logger.warn('No data reference established, experiment cannot be constructed')
                logger.warn(e)
                experiment = None
            if experiment is not None:
                TestConstruct_XML.make_happen(root, experiment, xml_file=xml_file, file_name=file_name)
        else:
            scell_file = open(file_name, 'r')
            text = scell_file.read()
            untagged_dict = untag(text)
            experiment, time_array = experiment_constructor(untagged_dict['Data_Reference'])
            untagged_dict['experiment'] = experiment
            TestConstruct.make_happen(untagged_dict)
            scell_file.close()
    # except Exception as e:
    #     logger.warn('Invalid file path')
    #     logger.warn(e)
    #     text = None
    #     scell_file = None


# for debugging purposes.  This module shouldn't do anything when run.
if __name__ == '__main__':
    read_from_file(r'C:\Users\Keeyan\Desktop\CCAM_Lab\sloppycell-git\source\temp\testFile.xml')
