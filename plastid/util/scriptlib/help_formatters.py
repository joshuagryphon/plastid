#!/usr/bin/env python
"""This module contains post-processors that reformats module docstrings for
use as command-line script help, by removing rich formatting, substitutions,
`reStructuredText`_ markup, and some `numpydoc`_ tokens.

See also
--------
:mod:`re`
    Python regular expressions module

`reStructuredText <http://docutils.sourceforge.net/rst.html>`_
    Specification for `reStructuredText`_, a markup language used by the
    `Sphinx`_ documentation engine, and used throughout the docstrings found
    in this package

`numpydoc <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
    A docstring standard and `Sphinx`_ extension that converts human-readable
    docstrings into valid `reStructuredText`_ 
"""
import re

def shorten_help(inp):
    """Pretty prints `numpydoc`_-formatted module, class, or function docstrings
    for use in command-line help, by removing some `reStructuredText`_ markup,
    and truncating the docstring at the appropriate `numpydoc`_ tokens.
    
    Parameters
    ----------
    inp : multi-line str
        Class or function docstring to format
    
    Returns
    -------
    str
        Cleaned helptext
    """
    inp = pyrst_pattern.sub(r"\g<spacing>\g<argument>",inp)
    inp = subst_pattern.sub(r"\g<1>",inp)
    inp = link_pattern.sub(r"\g<1>",inp)
    
    tokens = ["Parameters",
              "Returns",
              "Yields",
              "Raises",
              "Attributes"
              ]
    
    indices = [inp.find("%s" % X) for X in tokens]
    indices.extend([inp.find("    %s" % X.lower()) for X in tokens])
    
    for n,idx in enumerate(indices):
        if idx == -1:
            indices[n] = len(inp)
    
    return inp[:min(indices)].strip() + "\n"

def format_module_docstring(inp):
    """Pretty prints module docstrings for use in command-line help,
    by removing some `reStructuredText`_ markup, truncating the docstring
    at the first `numpydoc`_ token, and surrounding help text with separators.
    
    Parameters
    ----------
    inp : multi-line str
        Module docstring to format
    
    Returns
    -------
    str
        Formatted docstring
    """
    return _separator + "\n" + shorten_help(inp) + "\n" + _separator

pyrst_pattern = re.compile(r"(?P<spacing>^|\s+)(?::(?P<domain>[^:`<>]+))?:(?P<role>[^:`]*):`(?P<argument>[^`<>]+)(?: +<(?P<pointer>[^`]+)>)?`")
"""RegEx pattern that detects `reStructuredText`_ markup of python tokens
of the form ``:domain:role:`argument``` or simply ``:role:`argument```,
if the token is preceded by whitespace or begins a line.

Retrieves domain, role, and argument as named attributes if
:py:meth:`pyrst_pattern.groupdict` is called.
"""

subst_pattern = re.compile(r"\|([^|]*)\|")
"""RegEx pattern that matches `reStructuredText`_ substitution tokens
of form ``|substitution|``
"""

link_pattern = re.compile(r"`([^`<>]+)( <[^`]+>)?`_")
"""RegEx pattern that matches `reStructuredText`_ link references of forms ```Linkname`_``
and ```Link text <url>`_``
"""

_separator = "\n" + (78*"-") + "\n"
"""78-dash long text separator for separating help sections in command-line environments"""
