#!/usr/bin/env python
"""Miscellaneous functions that don't seem to fit elsewhere

:py:func:`guess_formatter`
    Convert string into primitive data type it probably encoded before ```str`` conversion

:py:func:`number`
    Convert string into numerical data type, prefering py:obj:`int` over py:obj:`float`
"""
import numpy

def guess_formatter(inp):
    """Guesses the format of input, preferring bools, ints, floats then strings
    Correctly parses nans and Infs. Converts None to nan


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
    Parses "nan", "Nan", "None", "none", "inf", and "-inf"


    Parameters
    ----------
    inp : str
        string input


    Returns
    -------
    float, numpy.nan, numpy.inf, or -numpy.inf, or str if no conversion found


    Raises
    ------
    ValueError : if type is non-numeric
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
