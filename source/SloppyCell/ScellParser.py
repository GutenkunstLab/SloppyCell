"""
Created 5/24/2017

Author @Keeyan

Parses .scell filetypes
"""
import logging
import re
logger = logging.getLogger('ScellParser')
logging.basicConfig()
# Defining global variables is bad practice so change this later @Keeyan
# TODO
tag_string = '<>'  # change things here to reformat Scell files
data_tag_string = '[]'


def tag(category, data):
    """
    Surrounds strings with the appropriate tags.
    Organizes data in a file to be read later.
    Returns two strings, a category and a data string.  Categories are contained 
    using the globally defined 'tag_string', and data are contained using the 
    'data_tag_string'.
    """
    print data
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

    # TODO: Use a keyedList to maintain consistency, also avoid directly modifying the dictionary object
    # Python is pass-by-object-reference, so modifying the object will work for references
    # elsewhere in the program, but reassignment will not register.
    # This could be a problem for maintenance/updates later, even if it works right now
    matching = [s for s in prelim_dict.keys() if "parameter" in s]
    for key in matching:
        params = prelim_dict[key]
        param_list = re.split(',', params)
        placeholder = {}
        for param in param_list:
            # This uses an exception to generate a parameter:value dictionary.  Probably not the best practice...
            # But it works.
            try:
                alpha = float(param)
                placeholder[beta] = alpha
            except ValueError:
                beta = param
        prelim_dict[key] = placeholder






def untag(file_string):
    """
    Turns scell file contents into a dictionary
    """

    # TODO: make regex function update dynamically with the globally defined tags (re.escape() could be helpful)
    return_dictionary = {}
    all_items = re.findall('<.*>\[.*\]', file_string)  # Makes a list of every item in the file
    for item in all_items:

        pair = re.split('>\[', item)

        return_dictionary[pair[0][1:]] = pair[1][:len(pair[1])-1]
    print return_dictionary
    retuple(return_dictionary)
    print return_dictionary




def write_to_file(package):
    target = open('testFile', 'w')
    for key in package.keys():
        a, b = tag(key, package[key])
        target.write(a + b + '\n')
    target.close()


def read_from_file(file):
    try:
        scell_file = open(file, 'r')
        text = scell_file.read()
        untag(text)
    except Exception as e:
        logger.warn('Invalid file directory')
        logger.warn(e)

if __name__ == '__main__':
    read_from_file(r'C:\Users\Keeyan\Desktop\CCAM_Lab\sloppycell-git\source\SloppyCell\testFile')
