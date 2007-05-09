import copy
import os

class KeyedList(list):
    def __init__(self, arg1 = None):
        # We stored both a dictionary mapping keys to indices and a list of
        #  keys for efficiency
        self.keyToIndex = {}
        self.storedKeys = []

        if hasattr(arg1, 'items'):
            items = arg1.items()
        else:
            items = arg1

        if items is not None:
            for (key, value) in items:
                self.set(key, value)

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def setOrder(self, order):
        if len(order) != len(self):
            raise ValueError, 'New order is of a different length!'

        oldOrder = copy.copy(self.keyToIndex)
        oldValues = copy.copy(self.values())
        self.storedKeys = []
        for (index, key) in enumerate(order):
            self.keyToIndex[key] = index
            self.storedKeys.append(key)
            self[index] = oldValues[oldOrder[key]]

    def reverse(self):
        self.setOrder(self.keys()[::-1])

    def __delitem__(self, index):
        index = index % len(self)

        list.__delitem__(self, index)
        del self.storedKeys[index]
        for key, ii in self.keyToIndex.items():
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
        return KeyedList(self.items())

    def __deepcopy__(self, memo):
        # XXX: Not handling recursion here
        return KeyedList(copy.deepcopy(self.items()))

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
            for key, value in other.items():
                self.set(key, value)
        else:
            if(len(self) != len(other)):
                raise ValueError, 'Other list not of same length!'
            for ii in range(len(self)):
                self[ii] = other[ii]

    def has_key(self, key):
        return self.keyToIndex.has_key(key)

    def setdefault(self, key, default=None):
        if not self.keyToIndex.has_key(key):
            self.set(key, default)

    def values(self):
        return self[:]

    def items(self):
        return zip(self.keys(), self.values())

    def __repr__(self):
        return 'KeyedList(%s)' % repr(self.items())

    def __str__(self):
        return os.linesep.join(['['] + 
                               [str(tup) + ',' for tup in self.items()] + [']'])

    #
    # list methods not supported
    #
    def  __add__(self, other):
        raise NotImplementedError
    def __iadd__(self, other):
        raise NotImplementedError
    def __imul__(self, factor):
        raise NotImplementedError
    def __mul__(self, factor):
        raise NotImplementedError
    def  __rmul__(self, factor):
        raise NotImplementedError
    def append(self, item):
        raise NotImplementedError
    def extend(self, other):
        raise NotImplementedError
    def insert(self, other):
        raise NotImplementedError
    def pop(self, other):
        raise NotImplementedError
    def remove(self, other):
        raise NotImplementedError
    def sort(self):
        raise NotImplementedError
