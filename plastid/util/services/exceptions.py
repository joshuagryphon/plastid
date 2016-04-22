#!/usr/bin/env python
"""This module contains custom exception and warning classes, implements 
a custom warning filter action, called `"onceperfamily"`, and monkey-patches
warning output to improve legibility.

Contents:

.. contents::
   :local:

The `onceperfamily` action
--------------------------
`onceperfamily` groups warning messages by families regular expressions,
and only prints the first warning instance that matches a given family's 
regular expression. In contrast, Python's native `once` action prints any string
literal once, even if it matches the same regex as another warning already given.

To use this action, use the following two functions:

  - :func:`filterwarnings` to create the warnings filter. Because :func:`filterwarnings`
    wraps Python's :func:`warnings.filterwarnings`, it may be used as a drop-in
    replacement for creation of any warnings filter.
    
  - :func:`warn` or :func:`warn_explicit`. Again, these are drop-in replacements
    for Python's :func:`warnings.warn` and :func:`warnings.warn_explicit` that
    additionally check the `onceperfamily` filters.
    
For convenience, there are also functions that both issue a warning, and create
a `onceperfamily` filter for it if one doesn't already exist:

  - :func:`warn_onceperfamily`
  - :func:`warn_explicit_onceperfamily`


Exception types
---------------
|MalformedFileError|
    Raised when a file cannot be parsed as expected, and
    execution must halt


Warning types
-------------
|ArgumentWarning|
    Warning for command-line arguments that:
    
      - are nonsenical, but recoverable
      - together might cause very slow execution
        (e.g. run would be optimized by other combinations)

|FileFormatWarning|
    Warning for slightly malformed but usable files

|DataWarning|
    Warning raised when:

      - data has unexpected attributes
      - data has nonsensical, but recoverable values for attributes
      - when values are out of the domain of a given operation,
        but skipping the operation or estimating the value is permissible


See also
--------
:mod:`warnings`
    Warnings module
"""
import re
import warnings
import inspect
import linecache
import textwrap
from plastid.util.io.filters import colored

_wrapper = textwrap.TextWrapper(break_long_words=False,width=77) 



#===============================================================================
# INDEX: Warning and exception classes
#===============================================================================

class MalformedFileError(Exception):
    """Exception class for when files cannot be parsed as they should be
    """
    
    def __init__(self,filename,message,line_num=None):
        """Create a |MalformedFileError|
        
        Parameters
        ----------
        filename : str
            Name of file causing problem
        
        message : str
            Message explaining how the file is malformed.
        
        line_num : int or None, optional
            Number of line causing problems
        """
        self.filename = filename
        self.msg      = message
        self.line_num = line_num
    
    def __str__(self):
        if self.line is None:
            return "Error opening file '%s': %s" % (self.filename, self.msg)
        else:
            return "Error opening file '%s' at line %s: %s" % (self.filename, self.line_num, self.msg)



class ArgumentWarning(Warning):
    """Warning for nonsensical but recoverable combinations of command-line arguments,
    or arguments that risk slow program execution"""
    pass


class FileFormatWarning(Warning):
    """Warning for slightly malformed but usable files"""
    pass


class DataWarning(Warning):
    """Warning for unexpected attributes of data.
    Raised when:

      - data has unexpected attributes
      - data has nonsensical, but recoverable values
      - values are out of the domain of a given operation, but execution
        can continue if the value is estimated or the operation skipped
    """



#===============================================================================
# INDEX: Plastid's extensions to Python warnings
#===============================================================================

pl_once_registry = {}
"""Registry of `onceperfamily` warnings that have been seen in the current execution context"""

pl_filters       = []
"""Plastid's own warnings filters, which allow additional actions compared to Python's"""

def filterwarnings(action,message="",category=Warning,module="",lineno=0,append=0):
    """Insert an entry into the warnings filter. Behaviors are as in :func:`warnings.filterwarnings`,
    except the additional action `'onceperfamily'` can be used to allow one warning per `family`
    of messages, specified by a regex. 
    
    This allows individual warnings to give more detailed information, without each being
    regarded as its own warning by Python's warning system (the defualt behavior of `'once'`).
    
    
    Parameters
    ----------
    action : str
        How the warning should be filtered. Accceptable values are "error",
        "ignore", "always", "default", 'module", "once", and "onceperfamily"
        
    message : str, optional
        str that can be compiled to a regex, used to detect warnings. If "onceperfamily"
        is chosen, only the first warning to give a string that matches the regex
        will be shown. For other actions, behaviors are as described in :mod:`warnings`.
        (Default: `""`, match any message)
        
    category : Warning or subclass, optional
        Type of warning. (Default: :class:`Warning`)

    module : str, optional
        str that can be compiled to a regex, limiting the warning behavior to modules
        that match that regex. (Default: `""`, match all modules)
        
    lineno : int, optional
        integer line used to specify warning in source code. If 0 (default), match
        all warnings regardless of line number.
        
    append : int, optional
        If 1, add warning to end of filter list. If 0 (default), insert warning at 
        beginning of filters list.
        
    
    See also
    --------
    warnings.filterwarnings
        Python's warnings filter
    """
    tup = (action,re.compile(message,re.I),category,re.compile(module),lineno)
    if action == "onceperfamily":
        if tup in pl_filters:
            return
        else:
            if append == 1:
                pl_filters.append(tup)
            else:
                pl_filters.insert(0,tup)
    else:
        warnings.filterwarnings(action,message=message,
                                category=category,module=module,
                                lineno=lineno,append=append)

def warn_onceperfamily(message,pattern=None,category=None,stacklevel=1):
    """Issue a warning and create a warning filter for that warning if it does not already exist
    
    Parameters
    ----------
    message : str
        Message of warning. Can be a string that compiles to a regex.
        Printed as warning text and used to create warning filter if `pattern`
        is `None`. 
    
    pattern : str or None, optional
        If not `None`, override message when generating warnings filter
        
    category: :class:`Warning`, or subclass, optional
        Type of warning
        
    stacklevel : int
        Ignored
        
    See also
    --------
    plastid.util.services.exceptions.filterwarnings
        plastid-specific warnings filters
    
    warnings.warn
        Python's warning system, which this wraps
    """
    if pattern is None:
        pattern = message
        filterwarnings("onceperfamily",message=pattern,category=category)
        
    warn(message,category=category,stacklevel=stacklevel)

def warn_explicit_onceperfamily(message,category,filename,lineno,pattern=None,
                                module=None,registry=None,module_globals=None):
    """Low-level interface to issue warnings, allowing `plastid`-specific warnings filters

    Parameters
    ----------
    message : str
        Message
    
    pattern : str or None, optional
        If not `None`, override message when generating warnings filter
        
    category: :class:`Warning`, or subclass, optional
        Type of warning

    filename : str
        Name of module from which warning is issued
    
    lineno : int, optional
        Line in module at which warning is called

    module : str, optional
        Module name
    
    registry : dict, optional
        Registry of ignore filters (see source code for :func:`warnings.warn_explicit`
    
    module_globals : dict, optional
        Dictionary of module-level variables
        
    
    See also
    --------
    plastid.util.services.exceptions.filterwarnings
        plastid-specific warnings filters

    warnings.warn_explicit
        Python's warning system, which this wraps
    """
    if pattern is None:
        pattern = message
        filterwarnings("onceperfamily",message=pattern,category=category)

    if module is None:# or module_globals is None:
        frame = inspect.currentframe()
        if frame is None:
            module = __name__
        else:
            try:
                module = inspect.getmodule(frame.f_back.f_code).__name__
            finally:
                del frame
        
    warn_explicit(pattern,category,filename,lineno,module=module,registry=registry,module_globals=module_globals)

def warn(message,category=None,stacklevel=1):
    """Issue a non-essential warning to users, allowing `plastid`-specific warnings filters
    
    Parameters
    ----------
    message : str
        Message
    
    category: :class:`Warning`, or subclass, optional
        Type of warning
        
    stacklevel : int
        Ignored
        
    See also
    --------
    plastid.util.services.exceptions.filterwarnings
        plastid-specific warnings filters
    
    warnings.warn
        Python's warning system, which this wraps
    """
    if category is None:
        category = UserWarning
        
    _, filename, lineno, _, _, _ = inspect.stack()[stacklevel]
    warn_explicit(message,category,filename,lineno,module=filename)

def warn_explicit(message,category,filename,lineno,module=None,registry=None,module_globals=None):
    """Low-level interface to issue warnings, allowing `plastid`-specific warnings filters

    Parameters
    ----------
    message : str
        Message
    
    category: :class:`Warning`, or subclass, optional
        Type of warning

    filename : str
        Name of module from which warning is issued
    
    lineno : int, optional
        Line in module at which warning is called

    module : str, optional
        Module name
    
    registry : dict, optional
        Registry of ignore filters (see source code for :func:`warnings.warn_explicit`
    
    module_globals : dict, optional
        Dictionary of module-level variables
        
    
    See also
    --------
    plastid.util.services.exceptions.filterwarnings
        plastid-specific warnings filters

    warnings.warn_explicit
        Python's warning system, which this wraps
    """
    global pl_once_registry
    if module is None: # or module_globals is None:
        frame = inspect.currentframe()
        if frame is None:
            module = __name__
        else:
            try:
                module = inspect.getmodule(frame.f_back.f_code).__name__
            finally:
                del frame

    for _, pat, filter_category, mod, filter_line in pl_filters:
        if pat.match(message) and issubclass(category,filter_category) and\
           (module is None or mod.match(module)) and\
           (filter_line == 0 or filter_line == lineno):
            
            tup =(pat.pattern,filter_category,mod,filter_line) 
            if tup in pl_once_registry:
                return
            else:
                pl_once_registry[tup] = 1
                break
            
    warnings.warn_explicit(message,category,filename,lineno,
                           module=module,registry=registry,
                           module_globals=module_globals)


def formatwarning(message,category,filename,lineno,file=None,line=None):
    """Wrapper to colorize warnings for readability. Overrides :func:`warnings.formatwarning`
    
    Parameters
    ----------
    message : str
        Warning message
        
    category : Warning
        Class (not instance) of warning
    
    filename : str
        Name of file calling warning
    
    lineno : int
        Line in file calling warning
    
    file : something implementing a `write` method 
        File to which output is written. Default: :obj:`sys.stderr`
        
    line : str
        Text of line in file calling warning. If `None`, `line` is taken
        to be line number `lineno` of `filename` 
    
    Returns
    -------
    str
        Pretty-printed warning message
    """
    sep     = colored("-"*75,color="cyan")
    message = str(message)
    if not "\n" in message:
        message = _wrapper.fill(str(message))
        
    message = colored(message,color="white",attrs=["bold"])
    name    = colored(category.__name__,color="cyan",attrs=["bold"])

    if line is None:
        numwidth = len(str(lineno+3))
        fmtstr   = "{0: >%ss} {1}" % (numwidth)
        lines    = []
        for x in range(max(0,lineno-2),lineno+3):
            tmpline = linecache.getline(filename,x).strip("\n")
            if tmpline:
                attrs = ["bold"] if x == lineno else []
                lines.append(fmtstr.format(colored(x,color="green",attrs=attrs),
                                           colored(tmpline,attrs=attrs)
                                           ))
        line = "\n".join(lines)

    filename = "in %s, line %s:" % (colored(filename,color="cyan"),lineno)

    ltmp = [sep,name,message,filename,"",line,"",sep,""]

    return "\n".join(ltmp)

    
warnings.formatwarning = formatwarning

