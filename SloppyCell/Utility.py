from __future__ import absolute_import
import pickle
import logging
import smtplib
import scipy
import numpy as np

from scipy import linspace, logspace

import random
import copy

def save(obj, filename):
    """
    Save an object to a file
    """
    f = open(filename, 'wb')
    pickle.dump(obj, f)
    f.close()

def load(filename):
    """
    Load an object from a file
    """
    f = open(filename, 'rb')
    obj = pickle.load(f)
    f.close()
    return obj

def eig(mat):
    """
    Return the sorted eigenvalues and eigenvectors of mat.
    """
    e, v = scipy.linalg.eig(mat)
    order = scipy.argsort(abs(e))[::-1]
    e = scipy.take(e, order)
    v = scipy.take(v, order, axis=1)

    return e, v

def bootstrap(data,num_iterates=100):
    """
    Use sampling with replication to calculate variance in estimates.
    """
    len_data = len(data)
    sampled_data = []
    for ii in range(num_iterates):
        sampled_data.append([random.choice(data) for ii in range(len_data)])

    return sampled_data
        
def enable_debugging_msgs(filename=None):
    """
    Enable output of debugging messages.

    If filename is None messages will be sent to stderr.
    """
    logging.root.setLevel(logging.DEBUG)

    if filename is not None:
        # Remove other handlers
        for h in logging.root.handlers:
            logging.root.removeHandler(h)

        # We need to add a file handler
        handler = logging.FileHandler(filename)
        # For some reason the default file handler format is different.
        # Let's change it back to the default
        formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
        handler.setFormatter(formatter)
        logging.root.addHandler(handler)
    logging.debug('Debug messages enabled.')

def disable_debugging_msgs():
    """
    Disable output of debugging messages.
    """
    logging.root.setLevel(logging.WARN)
    # Remove all other handlers
    for h in logging.root.handlers:
        logging.root.removeHandler(h)
    # Restore basic configuration
    logging.basicConfig()

def disable_warnings():
    np.seterr(over='ignore', divide='ignore', invalid='ignore',
                 under='ignore')
    logging.root.setLevel(logging.CRITICAL)

def enable_warnings():
    logging.root.setLevel(logging.WARNING)
    np.seterr(over='print', divide='print', invalid='print',
                 under='ignore')

class SloppyCellException(Exception):
    pass

class ConstraintViolatedException(SloppyCellException):
    """
    A ConstraintViolatedException is raised when a the value of a network
    constraint becomes False.
    """
    def __init__(self, time, constraint, message):
        self.constraint = constraint
        self.message = message
        self.time = time
    def __str__(self):
        return ('Violated constraint: %s at time %g. Additional info: %s.'
                %(self.constraint, self.time, self.message))

import SloppyCell.Redirector_mod
Redirector =  SloppyCell.Redirector_mod.Redirector

def combine_hessians(hesses, key_sets):
    """
    Combine a number of hessians (with possibly different dimensions and 
    orderings) into a single hessian.

    hesses    A sequence of hessians to combine
    key_sets  A sequence of sequences containing the parameter keys (in order)
              for each hessian.
    """

    # Get a list of all the possible parameter keys
    tot_keys = copy.copy(key_sets[0])
    keys_seen = set(tot_keys)
    for ks in key_sets[1:]:
        new_keys = [key for key in ks if key not in keys_seen]
        tot_keys.extend_keylist(new_keys)
        keys_seen|=set(new_keys)

    # Add all the appropriate hessian elements together
    key_to_index = dict(zip(tot_keys, range(len(tot_keys))))
    tot_hess = np.zeros((len(tot_keys), len(tot_keys)), scipy.float_)
    for hess, keys in zip(hesses, key_sets):
        for ii, id1 in enumerate(keys):
            tot_ii = key_to_index[id1]
            for jj, id2 in enumerate(keys):
                tot_jj = key_to_index[id2]
                tot_hess[tot_ii, tot_jj] += hess[ii, jj]

    return tot_hess, tot_keys
