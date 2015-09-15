#!/usr/bin/env python
"""Function decorators useful for scripts or analyses

Decorators
----------
:py:func:`catch_warnings`
    Catch warnings raised by wrapped function
    
:py:func:`deprecated`
    Function or class decorator. Wrapped functions or classes raise FutureWarnings
    when called or instantiated, respectively

:py:func:`parallelize`
    Parallelize the running of a single-parameter function over multiple instances
    of its parameter

:py:func:`catch_stdout`
    Redirect standard out from a wrapped function into a buffer

:py:func:`catch_stderr`
    Redirect standard error from a wrapped function into a buffer

:py:func:`in_separate_process`
    Run decorated function in a separate process, to force garbage collection
    of that function's memory contents when the function completes and reduce
    overall long-term memory usage.

:py:func:`notimplemented`
    Wrapped functions raise NotImplementedErrors. Use if committing incomplete code.

:py:func:`notused`
    No effects. Used for code annotation only

:py:func:`skip_if_abstract`
    Function decorator for unit tests. Wrapped methods will be skipped if
    they are called from a :py:class:`unittest.TestCase` with `'Abstract'`
    in its name, and run only in subclasses of the abstract :py:class:`unittest.TestCase`
    in which they are defined
"""
import functools
import os
import fcntl
import sys
import warnings
import types
import tempfile
from plastid.util.services.mini2to3 import get_func_code
                       
def notimplemented(func):
    """NotImplemented annotation decorator.
    Calls to functions annotated with this decorator raise
    |NotImplementedError|, which record attributes of function callers

    Parameters
    ----------
    func : function
        Function to wrap

    Returns
    -------
    function
        wrapped function
    """
    @functools.wraps(func)
    def new_func(*args,**kwargs):
        func_code = get_func_code(func)
        message = "NotImplementedException: call to unimplemented "+\
                  "function %s in module %s line %s" % (func.__name__,
                                                        func.__module__,
                                                        func_code.co_firstlineno + 1)
 
        raise NotImplementedError(message)
    
    return new_func
    
def notused(func):
    """Notused annotation decorator. Only for marking code.

    Parameters
    ----------
    func : function
        Function to wrap

    Returns
    -------
    function
        wrapped function
    """
    @functools.wraps(func)
    def new_func(*args,**kwargs):
        return func(*args,**kwargs)
    
    return new_func

def catch_warnings(simple_filter="ignore"):
    """Function factory producing function decorators that suppress warnings
    
    Parameters
    ----------
    simple_filter : str
        Warnings filter action. Quoted from :py:mod:`warnings`:
        
            ==============    ==================================================
            **Value**         **Disposition**
            --------------    --------------------------------------------------
            *error*           Turn warnings into exceptions
            *ignore*          Ignore all warnings
            *always*          Always print warnings
            *default*         Print first occurrence of each warning type,
                              for each location where warning is issued
            *module*          Print first occurrence of each warning type,
                              for each module where warning is issued
            *once*            Print first occurrence of matching warnings,
                              regardless of location
            ==============    ==================================================

        (source: :py:mod:`warnings`)

    
    Returns
    -------
    function
        Decorator function

        
    See also
    --------
    warnings
        Warnings module, especially sections on warnings filters
    """
    def decorator(func):
        """Function decorator that catches warnings of type %s
        
        Parameters
        ----------
        func : function
            Function to decorate
        
        Returns
        -------
        function
            Wrapped function
        """ % simple_filter
        @functools.wraps(func)
        def new_func(*args,**kwargs):
            with warnings.catch_warnings():
                warnings.simplefilter(simple_filter)
                result = func(*args,**kwargs)
            
            return result
        
        return new_func
    
    return decorator

def deprecated(func_or_class):
    """Deprecation annotation decorator for functions or classes. Wrapped
    functions or classes will raise FutureWarnings with called or instantiated,
    respectively.
    
    Parameters
    ----------
    func_or_class : function or class
        Function or class to deprecate

    Returns
    -------
    object
        wrapped function or class
    """
    # Based on useful hints from http://wiki.python.org/moin/PythonDecoratorLibrary
    if isinstance(func_or_class,types.FunctionType):
        @functools.wraps(func_or_class)
        def new_func(*args,**kwargs):
            message = "Call to deprecated function %s() which will be removed from future versions of module %s" % (func_or_class.__name__,
                                                                                                                    func_or_class.__module__)
            func_code = get_func_code(func_or_class)
            warnings.warn_explicit(message,
                                   category = DeprecationWarning,
                                   filename = func_code.co_filename,
                                   lineno = func_code.co_firstlineno + 1)
            return func_or_class(*args,**kwargs)
        
        return new_func
    # two tests here. Top line for Python 2.7. Bottom for 3.x
    elif (sys.version_info <= (3,) and (isinstance(func_or_class,types.ClassType) or isinstance(func_or_class,types.TypeType))) \
          or isinstance(func_or_class,type):
        old_init = func_or_class.__init__
        @functools.wraps(func_or_class.__init__)
        def new_func(self,*args,**kwargs):
            message = "Instantiation of deprecated class %s which will be removed from future versions of module %s" % (func_or_class.__name__,
                                                                                                                        func_or_class.__module__)
            warnings.warn_explicit(message,
                                   category = DeprecationWarning,
                                   filename = sys._getframe(1).f_globals["__name__"],
                                   lineno = sys._getframe(1).f_lineno )
            return old_init(self,*args,**kwargs)
        
        func_or_class.__init__ = new_func
        return func_or_class
    else:
        raise TypeError("Attempt to deprecate '%s', a %s. Only functions and classes can be deprecated." % (func_or_class.__name__,type(func_or_class)))

def skip_if_abstract(func):
    """Decorator function to keep :py:mod:`unittest` from running methods
    (defined in abstract classes) that are only intended to be run when
    inherited by fully-fleshed out subclasses. Wrapped methods are actually
    called from all non-abstract subclasses that inherit the method.
    
    Parameters
    ----------
    func : function
        Function that should only be run in a non-abstract subclass

    Returns
    -------
    function
        wrapped function
    """
    import unittest
    
    @functools.wraps(func)
    def new_func(*args,**kwargs):
        if "Abstract" in args[0].__class__.__name__:
            return unittest.skip("Skipping all tests from abstract class (don't worry, this is expected).")(func)
        else:
            return func(*args,**kwargs)
        
    return new_func

def catch_stderr(buf=None):
    """Function factory producing decorators that capture stderr to a buffer
    
    Parameters
    ----------
    buf : file
        Buffer that will hold captured stderr output. Must import
        ``write()`` and ``fileno()`` methods. If `None`, :py:obj:`os.devnull` will
        be used.
        
        
    Examples
    --------
    Create a ``pipe``, and use it to catch stderr::
    
        >>> import sys
        >>> import os
        >>>
        >>> def my_func():
        >>>    sys.stderr.write("some message")
        
        >>> read_end, write_end = os.pipe()
        >>> buf = os.fdopen(write_end,"w")

        >>> wrapped = catch_stderr(buf)(my_func)
        >>> wrapped() # generates stderr
    
        >>> buf.close()
        >>> captured = os.fdopen(read_end).read()
        >>> print(captured)
        # output from captured stderr here
        
        >>> captured.close() # remember to close!
    
    Returns
    -------
    function
        Function decorator
    """
    def decorator(func,buf=buf):
        """Decorator that suppresses standard error output from a function
        
        Parameters
        ----------
        func : function
            Function whose standard error will be captured
        
        Returns
        -------
        function
            wrapped function    
        """
        if buf is None:
            buf = open(os.devnull,"a")

        @functools.wraps(func)
        def new_func(*args,**kwargs):
            stderr_fd = os.dup(sys.__stderr__.fileno())
            new_fd  = buf.fileno()#tmpfile.fileno()
            os.dup2(new_fd,sys.stderr.fileno())
            try:
                result = func(*args,**kwargs)
            except BaseException as e:
                sys.stderr.flush() 
                os.dup2(stderr_fd,sys.stderr.fileno())
                raise(e)

            sys.stderr.flush() 
            os.dup2(stderr_fd,sys.stderr.fileno())
            return result

        return new_func

    return decorator

# cannot write unit tests for this because nose substitutes
# a StringIO object for sys.stdout, which messes up the file descriptor
# properties
def catch_stdout(buf=None):
    """Function factory producing decorators that capture stdout to a buffer
    
    Parameters
    ----------
    buf : file or None
        Buffer that will hold captured stdout output. Must import
        ``write()`` and ``fileno()`` methods. If `None`, :py:obj:`os.devnull` will
        be used.
        
        
    Examples
    --------
    Create a ``pipe``, and use it to catch stdout::
    
        >>> import sys
        >>> import os
        >>>
        >>> def my_func():
        >>>     sys.stdout.write("some message")
        
        >>> read_end, write_end = os.pipe()
        >>> buf = os.fdopen(write_end,"w")
        >>> 
        >>> wrapped = catch_stdout(buf)(my_func)
        >>> wrapped() # generates stdout
    
        >>> buf.close()
        >>> captured = os.fdopen(read_end).read()
        >>> print(captured)
        # output here
        
        >>> captured.close() # remember to close!
    
    Returns
    -------
    function
        Function decorator
    """
    def decorator(func,buf=buf):
        """Decorator that suppresses standard error output from a function
        
        Parameters
        ----------
        func : function
            Function whose standard error will be captured
        
        buf : file-like, optional
            Open stream, which **must** have a `fileno()` method. ``StringIO``
            objects will not work! (Default: `os.devnull`)
        
        Returns
        -------
        function
            wrapped function    
        """
        if buf is None:
            buf = open(os.devnull,"a")
            
        @functools.wraps(func)
        def new_func(*args,**kwargs):
            stdout_fd = os.dup(sys.stdout.fileno())
            new_fd  = buf.fileno()
            os.dup2(new_fd,sys.stdout.fileno())
            try:
                result = func(*args,**kwargs)
            except BaseException as e:
                sys.stdout.flush()
                os.dup2(stdout_fd,sys.stdout.fileno())
                raise(e)

            sys.stdout.flush()
            os.dup2(stdout_fd,sys.stdout.fileno())
            return result

        return new_func
    return decorator

def in_separate_process(func):
    """Decorator that runs a function in a separate process, to force garbage collection
    collection upon termination of that process, limiting long-term memory usage.
    
    Parameters
    ----------
    func : function
        Function to run in a separate process. Per multiprocessing spec,
        must be declared in global scope.
    
    Returns
    -------
    function
        wrapped function
    
    See Also
    --------
    multiprocessing
        Python :py:mod:`multiprocessing` module
    """
    @functools.wraps(func)
    def new_func(*args,**kwargs):
        import multiprocessing
        recv_pipe, send_pipe = multiprocessing.Pipe(False)
        
        def temp_func(conn,args,kwargs):
            result = func(*args,**kwargs)
            conn.send(result)
            conn.close()
        
        proc = multiprocessing.Process(target=temp_func,args=(send_pipe,args,kwargs))
        proc.start()
        proc.join()
        result = recv_pipe.recv()
        recv_pipe.close()
        
        return result
            
    return new_func

def parallelize(func):
    """Decorator to parallelize the running of a single-parameter function
    over multiple instances of its parameter(s)
     
    Parameters
    ----------
    func : function
        Function to parallelize. Must take one parameter. Per multiprocessing spec,
        must be declared in global scope.
     
    processes : int, optional
        Number of processes to use (Default: `4`)
     
    Returns
    -------
    function
        wrapped function
     
    See Also
    --------
    multiprocessing
        Python :py:mod:`multiprocessing` module
    """
    @functools.wraps(func)
    def new_func(args,processes=4,chunksize=None,**kwargs):
        import multiprocessing
        pf = functools.partial(func,**kwargs)
        pool = multiprocessing.Pool(processes=processes)
        pool_results = pool.map(pf,args,chunksize=chunksize)
        pool.close()
        pool.join()
         
        return pool_results
     
    message = """
     
    Notes
    -----
    #. This function has been parallelized with %s processes via
       :py:func:`plastid.util.services.decorators.parallelize`. Therefore,
       supply a lists of the function's parameter, rather than a single parameter.
    
    #. Output is not guaranteed to be sorted.
    
    #. This function additionally takes the keyword argument ``processes``,
       which determines how many processes it will use.
    """
     
    if new_func.__doc__ is not None:
        new_func.__doc__ += message
    else:
        new_func.__doc__ = message 
    return new_func

def skipdoc(func_or_class):
    """Instruct Sphinx not to generate documentation for ``func_or_class``

    Parameters
    ----------
    func_or_class : function or class
        Function or class not to document

    Returns
    -------
    object
        Wrapped function or class
    """
    func_or_class.plastid_skipdoc = True
    return func_or_class
