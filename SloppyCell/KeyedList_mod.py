"""
Holds the KeyedList class.
"""
from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import copy
import os

class KeyedList(list):
    def __init__(self, arg1 = None):
        # We stored both a dictionary mapping keys to indices and a list of
        #  keys for efficiency
        self.keyToIndex = {}
        self.storedKeys = []

        if hasattr(arg1, 'items'):
            items = list(arg1.items())
        else:
            items = arg1

        if items is not None:
            for (key, value) in items:
                self.set(key, value)

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def set_order(self, order):
        if len(order) != len(self):
            raise ValueError('New order is of a different length!')

        oldOrder = copy.copy(self.keyToIndex)
        oldValues = copy.copy(list(self.values()))
        self.storedKeys = []
        for (index, key) in enumerate(order):
            self.keyToIndex[key] = index
            self.storedKeys.append(key)
            self[index] = oldValues[oldOrder[key]]

    setOrder = set_order

    def reverse(self):
        self.set_order(list(self.keys())[::-1])

    def del_by_key(self, key):
        index = self.index_by_key(key)
        del self[index]

    def __delitem__(self, index):
        list.__delitem__(self, index)
        del self.storedKeys[index]
        for key, ii in list(self.keyToIndex.items()):
            if ii > index:
                self.keyToIndex[key] -= 1
            elif ii == index:
                del self.keyToIndex[key]

    def __delslice__(self, start, stop):
        # XXX: Not sure of behavior here
        if stop > len(self):
            stop = len(self)

        for index in range(start, stop)[::-1]:
            del self[index]

    def __copy__(self):
        instance = self.__new__(self.__class__)
        instance.__init__(list(self.items()))
        return instance

    def __deepcopy__(self, memo):
        # XXX: Not handling recursion here
        instance = self.__new__(self.__class__)
        instance.__init__(copy.deepcopy(list(self.items())))
        return instance

    #
    # Methods for manipulating by key.
    #
    def set(self, key, value):
        try:
            self[self.keyToIndex[key]] = value
        except KeyError:
            list.append(self, value)
            self.storedKeys.append(key)
            self.keyToIndex[key] = len(self)-1

    def index_by_key(self, key):
        return self.keyToIndex[key]

    def remove_by_key(self, key):
        del self[self.keyToIndex[key]]

    def get(self, key, default = None):
        try:
            return self[self.keyToIndex[key]]
        except KeyError:
            return default

    setByKey = set
    getByKey = get
    indexByKey = index_by_key
    removeByKey = remove_by_key

    def keys(self):
        return copy.copy(self.storedKeys)

    def update(self, other):
        if isinstance(other, dict) or isinstance(other, KeyedList):
            for key, value in list(other.items()):
                self.set(key, value)
        else:
            if(len(self) != len(other)):
                raise ValueError('Other list not of same length!')
            for ii in range(len(self)):
                self[ii] = other[ii]

    def has_key(self, key):
        return key in self.keyToIndex

    def setdefault(self, key, default=None):
        if key not in self.keyToIndex:
            self.set(key, default)

    def values(self):
        return self[:]

    def items(self):
        return list(zip(list(self.keys()), list(self.values())))

    def __repr__(self):
        return 'KeyedList(%s)' % repr(list(self.items()))

    def __str__(self):
        return os.linesep.join(['['] + 
                               [str(tup) + ',' for tup in list(self.items())] + [']'])

    #
    # list methods not supported
    #
    def __add__(self, other):
        new_kl = self.copy()
        new_kl += other
        return new_kl

    def __iadd__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError('Can only add another KeyedList to a KeyedList')
        for k,v in list(other.items()):
            if k not in self:
                self.set(k, v)
            else:
                raise ValueError('Addition would result in duplicated keys.')
        return self

    def extend_keylist(self, other):
        self += other

    def __imul__(self, factor):
        raise NotImplementedError
    def __mul__(self, factor):
        raise NotImplementedError
    def  __rmul__(self, factor):
        raise NotImplementedError

    def append(self, object):
        raise NotImplementedError

    def insert(self, index, object):
        raise NotImplementedError

    def insert_item(self, index, key, value):
        if key in self:
            raise ValueError('Insertion would result in duplicated key: %s.'
                             % str(key))
        list.insert(self, index, value)
        self.storedKeys.insert(index, key)
        for k, ii in list(self.keyToIndex.items()):
            if ii >= index:
                self.keyToIndex[k] += 1
        self.keyToIndex[key] = index

    def pop_value(self, index=-1):
        val = self[index]
        del self[index]
        return val
    pop = pop_value

    def pop_key(self, index=-1):
        key = self.storedKeys[index]
        del self[index]
        return key

    def pop_item(self, index=-1):
        item = (self.storedKeys[index], self[index])
        del self[index]
        return key

    def remove_by_value(self, other):
        index = self.index(other)
        del self[index]
    remove = remove_by_value

    def sort_by_key(self):
        """
        Sort based on key of each entry
        """
        keys = list(self.keys())
        keys.sort()
        self.set_order(keys)

    def sort_by_value(self):
        """
        Sort based on value of each entry
        """
        decorated = list(zip(list(self.values()), list(self.keys())))
        decorated.sort()
        sorted_keys = [k for (v,k) in decorated]
        self.set_order(sorted_keys)
    sort = sort_by_value
    
    def __contains__(self,key):
        return key in self.keys()
