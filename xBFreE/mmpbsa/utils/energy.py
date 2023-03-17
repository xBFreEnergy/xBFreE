# ##############################################################################
#                            GPLv3 LICENSE INFO                                #
#  Copyright (c) 2023.  Mario S. Valdes-Tresanco and Mario E. Valdes-Tresanco  #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify it     #
#  under the terms of the GNU General Public License version 3 as published    #
#  by the Free Software Foundation.                                            #
#                                                                              #
#  This program is distributed in the hope that it will be useful, but         #
#  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  #
#  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License    #
#  for more details.                                                           #
# ##############################################################################

import numpy as np
from math import sqrt
import pandas as pd
from typing import Union


class EnergyVector(np.ndarray):
    def __new__(cls, values=None, com_std=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        if isinstance(values, int):
            obj = np.zeros((values,)).view(cls)
        elif isinstance(values, (list, tuple, np.ndarray)):
            obj = np.array(values).view(cls)
        else:
            obj = np.array([]).view(cls)
        obj.com_std = com_std
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.com_std = getattr(obj, 'com_std', None)

    # This fix the pickle problem. Taken from
    # https://stackoverflow.com/questions/26598109/preserve-custom-attributes-when-pickling-subclass-of-numpy-array
    def __reduce__(self):
        # Get the parent's __reduce__ tuple
        pickled_state = super(EnergyVector, self).__reduce__()
        # Create our own tuple to pass to __setstate__
        new_state = pickled_state[2] + (self.__dict__,)
        # Return a tuple that replaces the parent's __setstate__ tuple with our own
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self.__dict__.update(state[-1])  # Update the internal dict from state
        # Call the parent's __setstate__ with the other tuple elements.
        super(EnergyVector, self).__setstate__(state[:-1])

    def stdev(self):
        return self.com_std or self.std()

    def sem(self):
        return float(self.std() / sqrt(len(self)))

    def semp(self):
        return float(self.stdev() / sqrt(len(self)))

    def append(self, values):
        return EnergyVector(np.append(self, values))

    def avg(self):
        return np.average(self)

    def corr_add(self, other):
        selfstd = self.com_std or float(self.std())
        comp_std = None
        if isinstance(other, EnergyVector):
            otherstd = other.com_std or float(other.std())
            comp_std = get_corrstd(selfstd, otherstd)
        return EnergyVector(np.add(self, other), comp_std)

    def corr_sub(self, other):
        self_std = self.com_std or float(np.asarray(self).std())
        comp_std = None
        if isinstance(other, EnergyVector):
            other_std = other.com_std or float(np.asarray(other).std())
            comp_std = get_corrstd(self_std, other_std)
        return EnergyVector(np.subtract(self, other), comp_std)

    def __add__(self, other):
        selfstd = self.com_std or float(self.std())
        comp_std = None
        if isinstance(other, EnergyVector):
            otherstd = other.com_std or float(other.std())
            comp_std = get_std(selfstd, otherstd)
        return EnergyVector(np.add(self, other), comp_std)

    def __sub__(self, other):
        self_std = self.com_std or float(np.asarray(self).std())
        comp_std = None
        if isinstance(other, EnergyVector):
            other_std = other.com_std or float(np.asarray(other).std())
            comp_std = get_std(self_std, other_std)
        return EnergyVector(np.subtract(self, other), comp_std)

    def __eq__(self, other):
        return np.all(np.equal(self, other))

    def __lt__(self, other):
        return np.all(np.less(self, other))

    def __le__(self, other):
        return np.all(np.less_equal(self, other))

    def __gt__(self, other):
        return np.all(np.greater(self, other))

    def __ge__(self, other):
        return np.all(np.greater_equal(self, other))

    def abs_gt(self, val):
        """ If any element's absolute value is greater than a # """
        return np.any(np.greater(np.abs(self), val))


def get_std(val1, val2):
    return sqrt(val1 ** 2 + val2 ** 2)


def get_corrstd(val1, val2):
    return sqrt(val1 ** 2 + val2 ** 2 - 2 * val1 * val2)


def calc_sum(vector1, vector2, mut=False) -> (float, float):
    """
    Calculate the mean and std of the two vector/numbers sum
    Args:
        vector1: EnergyVector or float
        vector2: EnergyVector or float
        mut: If mutant, the SD is the standard deviation of the array

    Returns:
        dmean: Mean of the sum
        dstd: Standard deviation
    """
    if isinstance(vector2, EnergyVector) and isinstance(vector1, EnergyVector):
        if mut:
            d = vector2 + vector1
            dmean = float(d.mean())
            dstd = float(d.std())
        else:
            dmean = float(vector2.mean() + vector1.mean())
            dstd = float(get_std(vector2.std(), vector1.std()))
    elif isinstance(vector2, EnergyVector) and isinstance(vector1, (int, float)):
        dmean = float(vector2.mean() + vector1)
        dstd = vector2.std()
    elif isinstance(vector2, (int, float)) and isinstance(vector1, EnergyVector):
        dmean = float(vector2 + vector1.mean())
        dstd = vector1.std()
    else:
        dmean = float(vector2 + vector1)
        dstd = 0.0
    return dmean, dstd


def multiindex2dict(p: Union[pd.MultiIndex, pd.Index, dict]) -> dict:
    """
    Converts a pandas Multiindex to a nested dict
    :parm p: As this is a recursive function, initially p is a pd.MultiIndex, but after the first iteration it takes
    the internal_dict value, so it becomes to a dictionary
    """
    internal_dict = {}
    end = False
    for x in p:
        # Since multi-indexes have a descending hierarchical structure, it is convenient to start from the last
        # element of each tuple. That is, we start by generating the lower level to the upper one. See the example
        if isinstance(p, pd.MultiIndex):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used. Only for 2 levels
            # pd.MultiIndex
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = [x[-1]]
            else:
                internal_dict[t].append(x[-1])
        elif isinstance(x, tuple):
            # This checks if the tuple x without the last element has len = 1. If so, the unique value of the
            # remaining tuple works as key in the new dict, otherwise the remaining tuple is used
            if len(x[:-1]) == 1:
                t = x[:-1][0]
                end = True
            else:
                t = x[:-1]
            if t not in internal_dict:
                internal_dict[t] = {x[-1]: p[x]}
            else:
                internal_dict[t][x[-1]] = p[x]
    if end:
        return internal_dict
    return multiindex2dict(internal_dict)


def flatten(dictionary, parent_key: list = False):
    """
    Turn a nested dictionary into a flattened dictionary
    :param dictionary: The dictionary to flatten
    :param parent_key:The accumulated list of keys
    :return: A flattened dictionary with key as tuples of nested keys
    """

    items = []
    for key, value in dictionary.items():
        new_key = parent_key + [key] if parent_key else [key]
        if isinstance(value, dict):
            items.extend(flatten(value, new_key).items())
        else:
            items.append((tuple(new_key), value))
    return dict(items)


def emapping(d):
    internal_dict = {}
    for k, v in d.items():
        if isinstance(v, dict):
            if v.values():
                if isinstance(list(v.values())[0], (dict, pd.DataFrame)):
                    internal_dict[k] = emapping(v)
                else:
                    internal_dict[k] = list(v.keys())
        elif isinstance(v, pd.DataFrame):
            internal_dict[k] = multiindex2dict(v.columns)
        else:
            internal_dict[k] = v
    return internal_dict
