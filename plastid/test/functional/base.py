#!/usr/bin/env python
"""Base classes for functional tests, which compare script output to reference output

These are implemented as :py:obj:`nose` test generators, so that tests may
be evaluated at the level of individual output files, as opposed to suites
of entire test runs. Thus, any test modules based upon this module will require
:py:obj:`nose` to run.

To use these tools, define a dictionary describing the test information
and the tests that need to be run. Then create a dummy function that calls
:py:func:`execute_helper`::

    #!/usr/bin/env python
    '''This module contains a functional test'''

    from module.I.want.to.test import main
    from plastid.test.functional.base import execute_helper
    from plastid.util.services.decorators import suppress_stderr
    
    my_test_info = { 
        "test_method"    : suppress_stderr(main), # hide output from screen
        "module_name"    : "module.I.want.to.test",
        "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
        "temp_file_path" : tempfile.mkdtemp(prefix="test_folder"),    
     }
     
    my_tests = [ ... ]
    
    def do_test():
        for x in execute_helper(my_test_info,my_tests):
            yield x
                
If run from the command-line, this script will generate modules for functional
tests of command-line scripts based upon a template resembling that above.
These modules may then be edited to flesh out the tests as necessary.
"""
import argparse
import sys
import shlex
import shutil
import inspect
import os
import nose

from nose.tools import assert_true
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import get_short_name, NullWriter
from plastid.bin.test_table_equality import main
from plastid.util.services.decorators import skip_if_abstract, catch_stderr
from pkg_resources import resource_filename, cleanup_resources

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))
table_test = catch_stderr()(main)

#===============================================================================
# INDEX: test generator functions
#===============================================================================

def check_one_file(test_info,ref_file,output_file,eq_args):
    """Check a single output file produced by a test run against a reference file 
    
    Parameters
    ----------
    test_info : dict
        Dictionary containing miscellaneous test information
    
    ref_file : str
        Name of reference file, excluding path, against which test output will be checked
    
    output_file : str
        Name of output file, excluding path
    
    eq_args : str
        Equality arguments determining how the reference and output file
        should be compared by :py:func:`table_test`
    """
    retval, failures = table_test([ref_file,output_file]+shlex.split(eq_args),verbose=True)
    msg = "\n    ".join(failures)
    err_message = "Unequal output in files %s vs %s, module %s:\n    %s" % (ref_file,output_file,test_info["module_name"],msg)
    nose.tools.assert_equal(retval,0,msg=err_message)

def build_test_lists(test_info,test):
    """Build a list of tests to be performed after one execution of a command-line script.
    
    Parameters
    ----------
    test_info : dict
        Dictionary containing miscellaneous test information
    
    test : list
        list of 3-element lists, each consisting of:
        
            1.  A reference filename against which test output will be compared
            
            2.  A filename of test output
            
            3.  A string of equality arguments determining how the reference
               and output file should be compared by :py:func:`table_test`
    
    Returns
    -------
    list
        list of:
        
        one 1-tuple describing arguments to past to the test method to run.
        This is always the first item in the list.
        
        The remaining items are one or more 5-tuples describing all the tests
        that must be run to cover each output file produced by the test method run.
        These contain:
        
            1.  the test function :py:func:`check_one_file`
            
            2.  the dictionary ``test_info``, describing file paths and other data
            
            3.  The file name of the reference output file for the test, 
                which the test output should match.
            
            4.  The corresponding test output filename for that test
            
            5.  Equality arguments for that test, to pass to :py:func:`table_test`
    """
    argstr, ref_files, test_files, eqarg_groups = test
    test_list = []
    test_list.append([argstr])
    for ref_file, output_file, eq_args in zip(ref_files,test_files,eqarg_groups):
        test_list.append((check_one_file, test_info, ref_file, output_file, eq_args))

    return test_list   

def execute_helper(test_info,tests):
    """Execute functional tests and evaluate all of their output files.
    This is a test generator function for use with :py:obj:`nose`. To useit,
    import it into a test module, and create a test function as follows::
    
        from plastid.test.functional.base import execute_helper
        my_test_info = { ... } # define this
        my_tests = [ ... ] # define tests. see Parameters section, below
        
        def do_test():
            for x in execute_helper(my_test_info,my_tests):
                yield x
    
    
    Parameters
    ----------
    test_info : dict
        Dictionary containing miscellaneous test information
    
    tests : list
        list of 4 element lists, each consisting of:
            
            1. Command-line arguments to pass to ``test_info["test_method"]``
            
            2. A list of reference files, which output files created by the test should match
             
            3. A list of output files created by the test
            
            4. A list of equality arguments (as strings) determining how each
               pair of reference and output files should be compared
               by :py:func:`table_test`
        
    Yields
    ------
    tuple
        5-tuples, collectively describing all the tests that must be run to cover
        every output file produced by every test run. Each tuple contains:
        
            1.  the test function :py:func:`check_one_file`
            
            2.  the dictionary ``test_info``, describing file paths and other data
            
            3.  The file name of the reference output file for the test, 
                which the test output should match.
            
            4.  The corresponding test output filename for that test
            
            5.  Equality arguments for that test, to pass to :py:func:`table_test`  
    """
    test_list = []
    for test in tests:
        _, ref_files, test_files, eqarg_groups = test
        assert len(eqarg_groups) == len(test_files) == len(ref_files)
        test_list.extend(build_test_lists(test_info,test))

    for test_item in test_list:
        # if len 1, item is an execute block, which will generate files to be tested next
        if len(test_item) == 1:
            # run test, but do it in try-catch so that if there is an execution error
            # subsequent tests will still run
            try:
                print("Executing test %s" % test_item[0])
                test_info["test_method"](shlex.split(test_item[0]))
            except Exception as e:
                # report that there was an error
                yield check_item_error, test_item,e 
        # otherwise, item is a test
        else:
            yield test_item

    if test_info["temp_file_path"] != "":
        shutil.rmtree(test_info["temp_file_path"])

def check_item_error(argstr,e):
    """Placeholder function used to raise exceptions when specific test executions fail
    
    Parameters
    ----------
    argstr : str
        String of arguments that were passed to test

    e : Exception
        Exception that was raised
    """
    assert_true(False,"Error executing test %s\n%s" % (argstr,e))

#===============================================================================
# INDEX: command-line program
#===============================================================================

def create_test_suite(module_names,base_path,overwrite=False,printer=NullWriter()):
    """Creates a :py:class:`unittest.TestCase` for the ``main()`` method of
    a command-line python module, and write these to modules
    
    Parameters
    ----------
    module_names : list
        List of module names, as strings, to create :py:class:`~unittest.TestCase` s for
    
    base_path : str
        Path where test modules will be written
    
    overwrite : bool, optional
        If True, existing files will be overwritten (Default: False)
    
    printer : file-like, optional
        Logger to which output will be written
    """
    for module_name in module_names:
        printer.write(module_name)
        short_name = get_short_name(module_name,separator="\.",terminator="")
        output_str = TEST_TEMPLATE.replace("${MODULE}",module_name)
        output_str = output_str.replace("${CAP_SHORT_NAME}",short_name.capitalize())
        output_str = output_str.replace("${SHORT_NAME}",short_name)
        output_fn  = os.path.join(base_path,"test_%s.py" % short_name)
        if overwrite or not os.path.exists(output_fn):
            printer.write("Writing %s" % output_fn)
            fout = open(output_fn,"w")
            fout.write(output_str)
            fout.close()

def main(argv=sys.argv[1:]):
    """Command-line program to generate functional test suites for command-line
    scripts from a template based upon the utilities in this module.
    
    For each test, output for given parameter sets is compared to reference output.
    
    Files will not be overwritten of they are already present.
    
    Parameters
    ----------
    argv : list
        list of command-line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("bin_folder",type=str,
                        help="Location of command-line modules to write tests for")
    parser.add_argument("output_folder",type=str,
                        help="Location of folder in which to write test modules")
    parser.add_argument("--overwrite",action="store_true",default=False,
                        help="If supplied, overwrite existing files")
    
    args = parser.parse_args(argv)
    modules   = [("plastid.bin.%s" % X).replace(".py","") \
                 for X in filter(lambda x: x.endswith(".py") and "__init__" not in x,
                                 os.listdir(args.bin_folder))]
    create_test_suite(modules,args.output_folder,overwrite=args.overwrite,printer=printer)
    
#===============================================================================
# INDEX : test template
#===============================================================================


TEST_TEMPLATE = '''#!/usr/bin/env python
"""Test suite for :py:mod:`${MODULE}`"""
import tempfile
import os

from nose.plugins.attrib import attr
from pkg_resources import resource_filename
from plastid.util.services.decorators import catch_stderr
from plastid.test.functional.base import execute_helper
from plastid.test.ref_files import REF_FILES, \\
                                       RPATH, \\
                                       REF_FILES, \\
                                       COUNT_OPTIONS, \\
                                       ANNOTATION_OPTIONS, \\
                                       MASK_OPTIONS
from ${MODULE} import main

test_info = {
    "test_method"    : catch_stderr()(main),
    "module_name"    : "${MODULE}",
    "ref_file_path"  : resource_filename("plastid","test/data/command_line"),
    "temp_file_path" : tempfile.mkdtemp(prefix="${SHORT_NAME}"),
}

# Define tests as tuples of:
# 1. Command-line style arguments to pass to ``main``
# 2. A list of reference files that output should be compared against
# 3. A list of output files created by running ``main`` with the arguments provided in (1)
# 4. A list of strings specifying how equality should be evaluated
${SHORT_NAME}_tests = [
    ("",[],[],[]),
    ("",[],[],[]),
    ("",[],[],[]),
    ("",[],[],[]),
    ("",[],[],[]),
]
"""Functional tests of :py:mod:`${MODULE}`.

Tests are specified as tuples of:

    1. Command-line style arguments to pass to :py:func:`main`

    2. A list of reference files that output should be compared against

    3. A list of output files created by running :py:func:`main`
       with the arguments provided in (1)

    4. A list of strings specifying how equality should be evaluated
"""


@attr(test="functional")
@attr(speed="slow")
def do_test():
    """Perform functional test for :py:mod`${MODULE}`"""
    for x in execute_helper(test_info,${SHORT_NAME}_tests):
        yield x
'''
"""Template for functional tests using test generators from :py:obj:`nose`"""
    
if __name__ == "__main__":
    main()
