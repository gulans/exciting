"""
Utilities for manipiulating dictionaries
"""
from typing import Callable, Iterator, Union
import json
import numpy as np
import copy
from collections.abc import MutableMapping, Mapping, Hashable
import sys


def common_iterable(obj: Union[dict, list]):
    """
    Create an iterable dict or tuple, such that one can iterate over a
    dictionary or list with the same syntax:

        for index in common_iterable(d):
            print('index = key or integer index:', index)
            print('d[index] = dictionary value or list element:', d[index])
    """
    if isinstance(obj, dict):
        return obj
    else:
        return (index for index, value in enumerate(obj))


def __container_converter(data: Union[list,dict]):
    """
    Mutates the input dictionary, converting string representations of numerical data into numerical data.

    :param dict data: Dictionary with string values.
    :return:dict new_data: Dictionary with all string literal values converted to numerical values.
    """
    np_convert = {np.float64: float,
                  np.int32: int}

    for element in common_iterable(data):
        if isinstance(data[element], dict) or isinstance(data[element], list):
           data[element] = __container_converter(data[element])
        elif isinstance(data[element], np.ndarray):
            data[element] = data[element].tolist()
        elif isinstance(data[element], (np.float64, np.int32)):
            value = data[element]
            data[element] = np_convert[type(value)](value)
        elif isinstance(data[element], Hashable):
            try:
                data[element] = json.loads(data[element])
            except:
                pass
        else:
            message = 'Type not converted by dict converter: ' + str(type(data[element]))
            sys.exit(message)

    return data


def container_converter(data: dict):
    """
    Mutates the input dictionary, converting string representations of numerical data into numerical data.

    :param dict data: Dictionary to be copied containing string values.
    :return: dict new_data: Copied dictionary with all string literal values converted to numerical values.
    """

    new_data = copy.deepcopy(data)
    return __container_converter(new_data)


def serialise_dict_values(d):
    """
    Mutate dictionary values that are objects, to dictionaries.
    Function works recursively, so no type-hinting.

    If any value within the dictionary is a set, it is converted to a list because
    common_iterable is not compatible with sets (sets cannot be indexed).

    If any value within the dictionary is a tuple, it is converted to a list because
    tuples are immutable, and therefore not consistent with an implementation that
    mutates the input.

    Sets and tuples are converted to lists regardless of whether they contain objects.

    :param dict input_dict: Dictionary.
    :param dict output_dict: Dictionary with any object values converted to dictionaries.
    """

    for index in common_iterable(d):
        # Can't iterate over a set with an index, and tuple cannot be mutated
        if isinstance(d[index], (set, tuple)):
            d[index] = list(d[index])

        # Class instances with the __dict__ attribute
        if '__dict__' in dir(d[index]):
            d[index] = vars(d[index])

        # Recursively apply routine to other containers
        elif not isinstance(d[index], Hashable):
            d[index] = serialise_dict_values(d[index])

        # pass over hashable values and np arrays
        else:
            pass

    return d


def get_hashable_entries(nested: dict) -> Iterator[tuple]:
    """
    Create a generator of all hashable values, removing nesting.

    Note, this will fail if a cross edge or a back edge is included in the dict
    i.e. a self-reference of the form:

        A = {key1: value, key2: A}

    However, for the use case of flattening parsed results, this should never be an issue.

    Reference: https://stackoverflow.com/questions/10756427/loop-through-all-nested-dictionary-values

    :param nested: Nested dictionary with a mix and hashable and non-hashable values.
    :return: Iterator[tuple] (key, value): Generator for tuples, containing keys and values.
    """
    for key, value in nested.items():
        if isinstance(value, abc.Mapping):
            yield from get_hashable_entries(value)
        else:
            yield key, value

