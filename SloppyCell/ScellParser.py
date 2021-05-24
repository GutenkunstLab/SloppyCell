"""
Created 6/14/2017

Author @Keeyan

Parses XML file, creates experiment
"""
from __future__ import print_function
from __future__ import absolute_import
import os
import csv
import logging
from xml.etree import ElementTree as ET
import pandas as pd
from . import TestConstruct_XML
from SloppyCell.ReactionNetworks import *

logger = logging.getLogger('ScellParser')
logging.basicConfig()




class GraphHandler:
    def __init__(self, graph_root):
        self.graph_root = graph_root


def experiment_constructor(data_file, sbml_reference):
    """
    Uses provided data file path to find appropriate csv file
    Extracts data and formats it to create an 'Experiment' object
    Returns the experiment object to be utilized in the main function
    """
    # TODO: Extend to allow multiple experiments to be constructed
    try:
        net = IO.from_SBML_file(sbml_reference)
        model_name = net.id
    except Exception as e:
        logger.warn("SBML reference not valid")
        logger.warn(e)
        return None
    expt = Experiment(os.path.splitext(os.path.basename(data_file))[0])  # Experiment named after data file
    # expt = Experiment('expt1')
    time_array = []
    expt_dict = {}
    if data_file.lower().endswith('.csv'):
        with open(data_file, 'rU') as csv_file:
            db = pd.read_csv(csv_file,skipinitialspace=True)

            new_col = db.columns.values
            new_col[0] = new_col[0].lower()
            new_col[1] = new_col[1].lower()
            db.columns =new_col
            ap = list(db['network'][db['network'].notnull()])
            db['network']=db['network'].fillna(method='ffill')
            for network in ap:
                return_dict =  combineSigma((db[db['network']==network]),network)
                expt_dict[network] = return_dict
        expt.set_data(expt_dict)
        try:
            time_array = list(db['time'])
        except Exception as e:
            time_array = list(db['Time'])
        return expt, time_array

    else:
        logger.warn('Selected data file type not supported.')

def combineSigma(df,network):
    lastcol = None
    return_dict = {}
    makelist = []
    for col in df:
        if 'sigma' in col.lower():
            dfa = df.dropna(subset=[col])
            dfa['mix']=list(zip(dfa[lastcol],dfa[col]))
            dfa_dict = pd.Series(dfa.mix.values, index=dfa.time).to_dict()
            if(len(dfa_dict) > 0):
                return_dict[lastcol] = dfa_dict

        lastcol = col

    return return_dict


def read_from_file(file_name, output_location=None):
    if output_location is not None:
        if not os.path.isdir(output_location):
            os.mkdir(output_location)

    if file_name.lower().endswith('.xml'):
        xml_file = ET.parse(file_name)
        root = xml_file.getroot()
        experiments = []
        try:
            sbml_reference = root.find("References").find("SBML").attrib['path']
            dir = os.path.dirname(os.path.abspath(file_name))
            #sbml_reference = os.path.relpath(sbml_reference,dir)
            sbml_reference = os.path.normpath(os.path.join(dir,sbml_reference))
            try:
                # TODO: Allow multiple data files for multiple experiments
                data_reference = root.find("References").findall('Data')
                dir = os.path.dirname(os.path.abspath(file_name))
                for ref in data_reference:
                    data_ref = ref.attrib['path']
                    #data_reference = os.path.relpath(data_reference, dir)
                    data_ref = os.path.normpath(os.path.join(dir, data_ref))
                    experiment, time_array = experiment_constructor(data_ref, sbml_reference)
                    experiments.append(experiment)
            except AttributeError as e:
                logger.warn('No data reference established, experiment cannot be constructed, SloppyCell will '
                            'have limited functionality.')
                logger.warn(e)
                experiment = None

        except TypeError as e:
            logger.warn('No sbml reference established, model cannot be made.')
            print(e)

        if experiments is not []:
            expt_dict = {}
            for experiment in experiments:
                expt_dict[experiment.GetName()] = experiment
            TestConstruct_XML.make_happen(root, experiment=expt_dict,
                                          xml_file=xml_file, file_name=file_name, sbml_reference=sbml_reference,
                                          output_location=output_location)
        else:
            TestConstruct_XML.make_happen(root, None, xml_file=xml_file, file_name=file_name,
                                         sbml_reference=sbml_reference, output_location=output_location)

# for debugging purposes.  This module shouldn't do anything when run.
if __name__ == '__main__':
    jak = r'C:\Users\ksg13004\Desktop\SloppyCell\sloppycell-git\Example\JAK-STAT\JAK-STAT.xml'
    tys = r'C:\Users\ksg13004\Desktop\SloppyCell\sloppycell-git\Example\Tyson_1991\Tyson1991.xml'
    pc12 = r'C:\Users\ksg13004\Desktop\SloppyCell\sloppycell-git\Example\PC12\PC12.xml'
    read_from_file(pc12)
