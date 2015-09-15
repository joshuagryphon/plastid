#!/usr/bin/env python
"""Tools for manipulating lists

Methods
-------
:py:func:`flatten_nested_lists_to_generator`
    Flatten nested lists, from left-to-right and top-to-bottom, into a generator

:py:func:`flatten_nested_lists_to_list`
    Flatten nested lists, from left-to-right and top-to-bottom, into a list

:py:func:`parse_list`
    Parse string representation of list into a list of parsed data (e.g.
    a list of numbers, a list of strings, et c)
"""

from plastid.util.services.misc import guess_formatter

def parse_list(inp):
    """Restore non-nested lists of Python primitives from string representation of list
    Additionally parses `numpy.nan`, `numpy.inf`, and `-numpy.inf` to correct values
    
    Parameters
    ----------
    inp : str
        String representation of list data

    Returns
    -------
    list

    Raises
    ------
    AssertionError
        if String does not represent a list

    Notes
    -----
    Complex types like nested lists or tuples will not be processed,
    and may not be parsed correctly.
    """
    assert inp[0] == '[' and inp[-1] == ']'
    return [guess_formatter(X.strip().strip("'")) for X in inp[1:-1].split(",")]

# Adapted from http://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists-in-python
def flatten_nested_lists_to_generator(l):
    """Flatten a tree of nested lists into a generator, from left-to-right and top-to-bottom.
        
    
    Parameters
    ----------
    l : list
        A tree as a list of lists
    
    Yields
    ------
    object
        Next item in flattened list
    """
    assert isinstance(l,list)
    for el in l:
        if isinstance(el, list):
            for sub in flatten_nested_lists_to_generator(el):
                yield sub
        else:
            yield el

def flatten_nested_lists_to_list(inp):
    """Flatten a tree of lists into a single list of items, from left-to-right and top-to-bottom.
    
    Parameters
    ----------
    inp : list
        A tree as a list of lists
    
    Returns
    -------
    list
        flattened list

    Notes
    -----
    Tuples and other sequences will not be flattened. They will remain in
    their native types
    """
    return list(flatten_nested_lists_to_generator(inp))
