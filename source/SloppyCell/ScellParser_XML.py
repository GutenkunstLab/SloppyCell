"""
Created 6/14/2017

Author @Keeyan

Parses XML file, creates experiment
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


def experiment_constructor(data_file, sbml_reference):
    """
    Uses provided data file path to find appropriate csv file
    Extracts data and formats it to create an 'Experiment' object
    Returns the experiment object to be utilized in the main function
    """
    try:
        net = IO.from_SBML_file(sbml_reference)
        model_name = net.id
    except Exception as e:
        logger.warn("SBML reference not valid")
        logger.warn(e)
        return None
    expt = Experiment(os.path.splitext(os.path.basename(data_file))[0])  # Experiment named after data file
    # expt = Experiment('expt1')
    if data_file.lower().endswith('.csv'):
        with open(data_file) as csv_file:
            reader_v = csv.DictReader(csv_file)
            field_list = reader_v.fieldnames
            # We assume the file is constructed with the format "Model,Time,Variable,Sigma,Variable,Sigma..."
            current_model = ''
            time = 0
            return_dictionary = {}
            switch = True
            data_name = []
            value = []
            sigma = []
            time_array = []
            for line in reader_v:
                for field in field_list:
                    decoded_field = field.decode("utf-8-sig").encode("utf-8")
                    if decoded_field.lower() == 'model':
                        if line[field] is not '':
                            current_model = line[field]
                    elif decoded_field.lower() == 'time':
                        time = int(line[field])
                    elif line[field] is not '':
                        if switch:
                            data_name.append(field)
                            value.append(line[field])
                            switch = False
                        elif not switch:
                            sigma.append(line[field])
                            switch = True
                for number in range(len(data_name)):
                    s = (time, data_name[number], float(value[number]), float(sigma[number]))
                    try:
                        return_dictionary[current_model][s[1]][time] = (s[2],s[3])
                        time_array.append(time)
                    except KeyError:
                        if current_model not in return_dictionary:
                            return_dictionary[current_model] = {}
                        if s[1] not in return_dictionary[current_model]:
                            return_dictionary[current_model][s[1]] = {}
                        return_dictionary[current_model][s[1]][time] = (s[2], s[3])
                        time_array.append(time)
                data_name = []
                value = []
                sigma = []
            expt.set_data(return_dictionary)
            return expt, time_array

    else:
        logger.warn('Selected data file type not supported.')


def read_from_file(file_name):
    if file_name.lower().endswith('.xml'):
        xml_file = ET.parse(file_name)
        root = xml_file.getroot()
        experiment = None
        try:
            sbml_reference = root.find("References").find("SBML").attrib['path']
            try:
                data_reference = root.find("References").find("Data").attrib["path"]
                experiment, time_array = experiment_constructor(data_reference, sbml_reference)
            except AttributeError as e:
                logger.warn('No data reference established, experiment cannot be constructed')
                logger.warn(e)
                experiment = None

        except TypeError as e:
            logger.warn('No sbml reference established, model cannot be made.')
            print e
        if experiment is not None:
            TestConstruct_XML.make_happen(root, experiment, xml_file=xml_file, file_name=file_name,
                                          sbml_reference=sbml_reference)

# for debugging purposes.  This module shouldn't do anything when run.
if __name__ == '__main__':
    read_from_file(r'C:\Users\Keeyan\Desktop\CCAM_Lab\sloppycell-git\source\temp\jak_stat_example.xml')