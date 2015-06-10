#!/usr/bin/env python
"""This module contains post-processors that reformats module docstrings for
use as command-line script help, by removing rich formatting, substitutions,
and reStructuredText tokens.
"""
import re

def shorten_help(inp):
    """Pretty prints numpydoc-formatted docstrings for use in command-line help,
    by removing some reStructuredText markup, and truncating the docstring
    at the appropriate numpydoc tokens.
    
    Parameters
    ----------
    inp : str
        Class or function docstring to format
    
    Returns
    -------
    str
        Cleaned helptext
    """
    inp = pyrst_pattern.sub(r"\g<spacing>\g<argument>",inp)
    inp = subst_pattern.sub(r"\g<1>",inp)
    
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
    by removing some reStructuredText markup, truncating the docstring
    at the first numpydoc token, and surrounding help text with separators.
    
    Parameters
    ----------
    inp : str
        Module docstring to format
    
    Returns
    -------
    str
        Formatted docstring
    """
    return _separator + shorten_help(inp) + _separator

#TODO: add support for :domain:role:`Text <argument>`
#TODO: add support for links `Link`_
pyrst_pattern = re.compile(r"(?P<spacing>^|\s+)(?::(?P<domain>py))?:(?P<role>[^:]*):`(?P<argument>[^`]*)`")
"""RegEx pattern that detects reStructuredText markup of python tokens
of the form ``:domain:role:`argument``` or simply ``:role:`argument```,
if the token is preceded by whitespace or begins a line.

Retrieves domain, role, and argument as named attributes if
:py:meth:`pyrst_pattern.groupdict` is called.

See Also
--------
re
    Python regular expressions module
"""

subst_pattern = re.compile(r"\|([^|]*)\|")
"""RegEx pattern that matches reStructuredText substitution tokens

See Also
--------
re
    Python regular expressions module
"""

_separator = "\n" + (78*"-") + "\n"
"""78-dash long text separator for separating help sections in command-line environments"""
