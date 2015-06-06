#!/usr/bin/env python
"""Utilities for converting colors between formats

:py:func:`get_rgb255_from_str`
    Converts input of RGB hex strings of type #RRGGBB to rgb tuples

:py:func:`get_str_from_rgb255`
    Converts RGB tuples to RGB hex strings of type #RRGGBB
"""
def get_rgb255_from_str(inp):
    """Converts input of RGB hex strings of type #RRGGBB to rgb tuples

    Parameters
    ----------
    inp : str
        RGB hex string of form '#NNNNNN' where 0 <= 0xNN <= 255


    Returns
    -------
    tuple<int 0..255,int 0..255,int 0..255>


    Raises
    ------
    ValueError if values are out of range
    """
    r = int("0x%s" % inp[1:3],16)
    g = int("0x%s" % inp[3:5],16)
    b = int("0x%s" % inp[5:],16)
    if len(inp) > 7 or r > 255 or g > 255 or b > 255:
        raise ValueError("Cannot parse non-hex string (should be '#NNNNNN' got '%s')" % inp)
    return (r,g,b)

def get_str_from_rgb255(inp):
    """Converts RGB tuples to RGB hex strings of type #RRGGBB
    
    Parameters
    ----------
    input : tuple<int 0..255,int 0..255,int 0..255>


    Returns
    -------
    str
        RGB hex string of form '#NNNNNN'
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
