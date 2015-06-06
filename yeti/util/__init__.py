#!/usr/bin/env python
"""Miscellaneous utilities

    ========================================   ======================================================================================================================
    Sub-package                                Contents
    ----------------------------------------   ----------------------------------------------------------------------------------------------------------------------
    :py:obj:`~yeti.util.io`             Wrappers for various file I/O operations
    :py:obj:`~yeti.util.services`       Functions that operate on miscellaneous basic data types (e.g. sets, lists, et c), function decorators, and exceptions
    ========================================   ======================================================================================================================


    ========================================   ======================================================================================================================
    Module                                     Contents
    ----------------------------------------   ----------------------------------------------------------------------------------------------------------------------
    :py:mod:`~yeti.util.array_table`     DataFrame-like object (deprecated)
    :py:mod:`~yeti.util.unique_fifo`    FIFO that maintains a unique collection of members. If a member already present in the FIFO as added to the FIFO,
                                               it is moved from its original position to the end FIFO, instead of retaining its original position and being added
                                               again at the end.
    ========================================   ======================================================================================================================



"""
