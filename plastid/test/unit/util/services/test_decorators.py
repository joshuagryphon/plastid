#!/usr/bin/env python
import sys
import os
import fcntl
import types
import warnings
from nose.plugins.attrib import attr
from nose.tools import assert_equal, assert_true, assert_greater_equal, assert_raises, assert_set_equal, assert_not_equal
from plastid.util.services.decorators import notimplemented, \
                                                 deprecated, \
                                                 parallelize, \
                                                 in_separate_process, \
                                                 catch_stderr, \
                                                 catch_stdout, \
                                                 catch_warnings

#===============================================================================
# INDEX: functions and classes that will be decorated
#===============================================================================

def stderr_func(msg):
    sys.stderr.write(msg)
    return msg

def stdout_func(msg):
    sys.stdout.write(msg)
    return msg

def util_func(x):
    """Square numbers and return process in which function was run
    
    Parameters
    ----------
    x : int or float

    Returns
    -------
    int or float
        Squared value of ``x``
    
    int
        Process ID in which function was run
    """
    return x**2, os.getpid()

def func_that_warns():
    warnings.warn("Some warning",UserWarning)

def get_pipes():
    readfd, writefd = os.pipe()
    fcntl.fcntl(readfd,fcntl.F_SETFL,os.O_NONBLOCK)
    fcntl.fcntl(writefd,fcntl.F_SETFL,os.O_NONBLOCK)
    return os.fdopen(readfd,"r"), os.fdopen(writefd,"w")

class UtilClass(object):
    def __init__(self,tmp):
        self.name = tmp
    
    def get_name(self):
        return self.name

#===============================================================================
# INDEX: unit tests
#===============================================================================

# stdout/err redirection -------------------------------------------------------
STDERR_FD = sys.stderr.fileno()

if sys.version_info[0] == 2:
    EMPTY_BUFFER_ERROR = IOError
else:
    # Python 3.x codecs returns None from empty non-blocking buffers
    # which causes a TypeError to be raised instead of 
    # the IOError raised in Python 2.x
    EMPTY_BUFFER_ERROR = TypeError

@attr(test="unit")
def test_catch_stderr_doesnt_print_without_buffer():
    # spy on `inner` by making sure there is nothing written to stderr
    outer_reader, outer_writer = get_pipes()
    message = "this is a test"
    sys.stderr = sys.__stderr__
    
    @catch_stderr(outer_writer)
    def inner():
        # make sure value is returned from wrapped function
            wrapped = catch_stderr()(stderr_func)
            msg = wrapped(message)
            assert_equal(msg,message)
    
    inner()
    outer_writer.flush()
    outer_writer.close()

    # make sure no message made it out of `inner`
    # with nothing to read in pipe, outer_reader should raise IOError
    assert_raises(EMPTY_BUFFER_ERROR,outer_reader.read)
    outer_reader.close()
    
@attr(test="unit")
def test_catch_stderr_doesnt_print_with_buffer_but_catches_in_buffer():
    # spy on `inner` by making sure there is nothing written to stderr
    # but make sure message is found in inner readre
    outer_reader, outer_writer = get_pipes()
    message = "this is a test"
    sys.stderr = sys.__stderr__

    @catch_stderr(outer_writer)
    def inner():
        inner_reader, inner_writer = get_pipes()
        wrapped = catch_stderr(inner_writer)(stderr_func)
        
        # make sure value is returned from wrapped function
        msg = wrapped(message)
        assert_equal(message,msg)
        
        # make sure we caught entire message from stderr
        inner_writer.flush()
        inner_writer.close()
        assert_equal(message,inner_reader.read())
        inner_reader.close()
    
    inner()
    outer_writer.flush()
    outer_writer.close()

    # make sure no message made it out of `inner`
    # with nothing to read in pipe, outer_reader should raise IOError
    assert_raises(EMPTY_BUFFER_ERROR,outer_reader.read)
    outer_reader.close()



# catch warnings    ------------------------------------------------------------

@attr(test="unit")
def test_catch_warnings_catches_warnings():
    ign_func = catch_warnings("ignore")(func_that_warns)

    num = 5
    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        for x in range(num):
            func_that_warns()
    
    # make sure warning is issued with vanilla function
    assert_equal(len(warns),num)

    num = 5
    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        for x in range(num):
            ign_func()
    
    # make sure warning is caught with wrapped function
    assert_equal(len(warns),0)


# notimplemented    ------------------------------------------------------------

@attr(test="unit")
def test_notimplemented_raises_exception():
    my_func = notimplemented(util_func)
    assert_true(isinstance(my_func,types.FunctionType))
    
    assert_raises(NotImplementedError,my_func,5)

# deprecated    ----------------------------------------------------------------

@attr(test="unit")
def test_deprecated_function_raises_warning_only_once():
    num = 5
    my_func = deprecated(util_func)
    assert_true(isinstance(my_func,types.FunctionType))
    
    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        for x in range(num):
            assert_equal(my_func(x),util_func(x))
            my_func(x)
    
    # make sure warning is issued only once
    assert_equal(len(warns),1)

@attr(test="unit")
def test_deprecated_class_raises_warning():
    reg_obj = UtilClass("my_object")
    dep_class = deprecated(UtilClass)

    with warnings.catch_warnings(record=True) as warns:
        warnings.simplefilter("always")
        dep_obj = dep_class("my_object")
        assert_true(isinstance(dep_obj,UtilClass))
    
    # make sure warning is issued
    assert_equal(len(warns),1)
    
    # make sure wrapped class behaves as it should
    assert_equal(reg_obj.get_name(),dep_obj.get_name())


# parallelize or in other processes --------------------------------------------

@attr(test="unit")
def test_parallelize_spawns_processes_and_gets_correct_ansswer():
    x = range(500)
    my_func = parallelize(util_func)
    outer_vals, outer_pids = zip(*[util_func(X) for X in x])
    inner_vals, inner_pids = zip(*my_func(x))

    assert_set_equal(set(outer_vals),set(inner_vals))
    assert_not_equal(set(outer_pids),set(inner_pids))

@attr(test="unit")
def test_in_separate_process_spawns_process_and_gets_correct_ansswers():
    my_func = in_separate_process(util_func)
    assert_true(isinstance(my_func,types.FunctionType))
    
    for x in range(100):
        res_out, outer_pid = my_func(x)
        res_in,  inner_pid = util_func(x)        
        assert_equal(res_out,res_in)
        assert_not_equal(outer_pid,inner_pid)
        
