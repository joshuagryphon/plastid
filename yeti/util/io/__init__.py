#!/usr/bin/env python
"""Utilities for handling and wrapping various I/O operations.

Package overview
================

    =======================================  ======================================
    **Package module**                       **Contents**
    ---------------------------------------  --------------------------------------
    :py:mod:`~plastid.util.io.binary`           Tools for unpacking binary values
                                             into named dictionaries 
    :py:mod:`~plastid.util.io.filters`          Wrappers that transform input or
                                             output during file I/O (e.g. remove
                                             comments, skip blank lines, perform
                                             arbitrary function before reading
                                             or writing).
    :py:mod:`~plastid.util.io.openers`          Wrappers for opening and closing files
    =======================================  ======================================

"""

