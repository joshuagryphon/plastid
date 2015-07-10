#!/usr/bin/env python
"""Miscellaneous functions that don't seem to fit elsewhere

:py:func:`guess_formatter`
    Convert string into primitive data type it probably encoded before `str`
    was called on it

:py:func:`number`
    Convert string into numerical data type, trying `bool`, `int`, `float`,
    then `str`
"""
import numpy

def guess_formatter(inp):
    """Guesses the format of input, trying `bool`, `int`, `float`, then `str`.
    Correctly parses `nan`s and `Inf`s. Converts `None` to `nan`


    Parameters
    ----------
    inp : str
        input


    Returns
    -------
    boolean, number, or string
    """
    if inp.lower() == "true":
        return True
    elif inp.lower() == "false":
        return False
    else:
        try:
            return number(inp)
        except ValueError:
            return str(inp)
    
def number(inp):
    """Parses numbers from strings, preferring int over float.
    Parses `nan`, `Nan`, `None`, `none`, `inf`, and `-inf`


    Parameters
    ----------
    inp : str
        string input


    Returns
    -------
    float, numpy.nan, numpy.inf, or -numpy.inf, or str if no conversion found


    Raises
    ------
    ValueError
        if `inp` cannot be converted to a number
    """
    if inp in ("nan","NaN","na","None","none"):
        return numpy.nan
    elif inp in ("inf","Inf"):
        return numpy.inf
    elif inp in ("-inf","-Inf"):
        return -numpy.inf
    else:
        try:
            # note: in python bools are also ints! 
            # isinstance(True,int) == True
            val = int(inp)
        except ValueError:
            val = float(inp)
        return val
