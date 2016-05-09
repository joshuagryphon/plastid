#!/usr/bin/env python
"""Utilities for manipulating colors and converting between representations
as RGB hex strings, tuples of floats from 0.0 to 1.0, and tuples of ints
from 0.0 to 255.
"""
import numpy
from matplotlib.colors import colorConverter

process_black = "#222222"

def get_rgb255(inp):
    """Fetch `r,g,b` values where `r,g,b` are integers ranging from 0 to 255

    Parameters
    ----------
    inp : *RGB* or *RGBA* sequence or str
        Can be color hex string, like "#007ADF", a matplotlib color letter like "r",
        a matplotlib color name like "black, et c."
        See :meth:`matplotlib.colors.colorConverter.to_rgb`

    Returns
    -------
    :class:`numpy.ndarray`
        Numpy array of `r,g,b` tuples where `r,g,b` take integer values from 0 to 255
    """
    data = colorConverter.to_rgba_array(inp)
    data = (255*data[:,:3]).round().astype(int).ravel()

    return data

def get_str_from_rgb(inp):
    """Converts RGB tuples of floats from between 0.0 and 1.0 to RGB hex strings of type #RRGGBB
    
    Parameters
    ----------
    input : tuple
        Tuple of r,g,b values in range from 0.0 to 1.0

    Returns
    -------
    str
        RGB hex string of form '#NNNNNN'

    Raises
    ------
    ValueError if values are out of range
    """    
    return get_str_from_rgb255((255*numpy.array(inp).round()).astype(int))
    
def get_str_from_rgb255(inp):
    """Converts RGB tuples of ints from between 0 and 255 to RGB hex strings of type `#RRGGBB`
    
    Parameters
    ----------
    input : tuple
        Tuple of r,g,b values in range from 0 to 255


    Returns
    -------
    str
        RGB hex string of form '#NNNNNN'
        
    Raises
    ------
    ValueError if values are out of range        
    """
    err_msg = "Cannot convert malformed tuple to string. Should be (0..255,0..255,0..255). Is: (%s)" % ",".join([str(X) for X in inp])
    try:
        assert len(inp) == 3
        for x in inp:
            assert x >= 0 and x <= 255
        stmp = "#{:0>2X}{:0>2X}{:0>2X}".format(*inp)
        assert len(stmp) == 7
        return stmp
    except IndexError:
        raise ValueError(err_msg)
    except AssertionError:
        raise ValueError(err_msg)
    except TypeError:
        raise ValueError(err_msg)

def lighten(data,amt=0.10,is255=False):
    """Lighten a vector of colors by fraction `amt` of remaining possible intensity.
    
    New colors are calculated as::
     
        >>> new_colors = data + amt*(1.0-data)
        >>> new_colors[:,-1] = 1 # keep all alpha at 1.0
    
    Parameters
    ----------
    data : matplotlib colorspec or sequence of colorspecs
        input color(s)
    
    amt : float, optional
        Percentage by which to lighten `r`, `g`, and `b`. `a` remains unchanged
        (Default: 0.10)

    is255 : bool, optional
        If `True`, rgb values in `data` are assumed to be tween 0 and 255 
        rather than 0.0 and 1.0. In this case, return values will also 
        be between 0 and 255.
        
    Returns
    -------
    numpy.ndarray
        Lightened version of data
    """
    data = colorConverter.to_rgba_array(data)

    new_colors = data + amt*(1.0-data)
    if is255:
        new_colors = (255*new_colors).round()
    
    new_colors[:,-1] = data[:,-1]

    return new_colors

def darken(data,amt=0.10,is255=False):
    """Darken a vector of colors by fraction `amt` of current intensity.
    
    Parameters
    ----------
    data : matplotlib colorspec or sequence of colorspecs
        input color(s)

    amt : float, optional
        Percentage by which to darken `r`, `g`, and `b`. `a` remains unchanged
        (Default: 0.10)
        
    Returns
    -------
    numpy.ndarray
        Lightened version of data
    """
    data = colorConverter.to_rgba_array(data)

    new_colors = (1.0-amt)*data
    if is255:
        new_colors = (255*new_colors).round()

    new_colors[:,-1] = data[:,-1]

    return new_colors

