#!/usr/bin/env python
"""This module provides a thin compatibility layer between Python 2.7 and 3.x.
Various objects are aliased as follows:

    ===================    ===========================   =======================
    **Exported object**     **Points to in 2.x**         **Points to 3.x**
    -------------------    ---------------------------   -----------------------
    ``cStringIO``          :mod:`StringIO`               :mod:`io`
    ``StringIO``           :mod:`cStringIO`              :mod:`io`
    ``xrange``             :func:`xrange`                :func:`range`
    ``quote``              :func:`urllib.quote`          :func:`urllib.parse.quote`
    ``unquote``            :func:`urllib.quote`          :func:`urllib.parse.unquote`
    ``quote_plus``         :func:`urllib.quote_plus`     :func:`urllib.parse.quote_plus`
    ``unquote_plus``       :func:`urllib.unquote_plus`   :func:`urllib.parse.unquote_plus`
    ``ifilter``            :func:`itertools.ifilter`     :func:`filter`
    ===================    ===========================   =======================


Also, one function is defined:

    ========================    ================================================
    **Function**                **Action**
    ------------------------    ------------------------------------------------
    :func:`get_func_code`       Retrieves ``function.func_code`` in 2.x, ``function.__code__`` in 3.x
    ========================    ================================================

"""
import sys
import urllib
import itertools

if sys.version_info >= (3,):
    import io as cStringIO
    import io as StringIO
    xrange = range
    quote   = urllib.parse.quote
    unquote = urllib.parse.unquote
    quote_plus   = urllib.parse.quote_plus
    unquote_plus = urllib.parse.unquote_plus
    ifilter = filter
    
    # function code
    _func_code_attr = "__code__"
else:
    import cStringIO
    import StringIO
    xrange = xrange
    quote   = urllib.quote
    unquote = urllib.unquote
    quote_plus   = urllib.quote_plus
    unquote_plus = urllib.unquote_plus
    ifilter = itertools.ifilter
    
    # function code
    _func_code_attr = "func_code"


def get_func_code(func):
    """Return function code attribute for function ``func``
    
        ==========     ===================
        Version        Function attribute
        ----------     -------------------
        2.x            ``func_code``
        3.x             ``_code``
        ==========     ===================


    Parameters
    ----------
    func : function
        Query function
    
    Returns
    -------
    code
        Function code
    """
    return getattr(func,_func_code_attr)