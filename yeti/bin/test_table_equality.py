#!/usr/bin/env python
"""Regression-testing script designed to test equality between a newly-generated
file and a reference file that is intended to contain the same data. Rows and 
columns are not expected to be in the same order. Float values are only required
to be equal within a user-specified tolerance. `NaN` values evaluate as equal
if and only if they occur in the same cell after sorting rows and columns.
Finally, specific columns may be excluded by name (or number, if there is no
header row).

Exit status is 0 if files are identical, 1 otherwise 
"""
from yeti.util.array_table import ArrayTable, equal_enough
from yeti.util.io.openers import opener, get_short_name, NullWriter
from yeti.util.io.filters import NameDateWriter, CommentReader
from pandas.util.testing import assert_frame_equal
import sys
import pandas as pd
import argparse
import inspect

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

NUMERIC_DTYPES = "biufc"

def test_dataframe_equality(df1,df2,tol=1e-8,sort_columns=[],printer=NullWriter(),print_verbose=False,return_verbose=False):
    """Test equality of dataframes over multiple columns, with verbose output.
    If `NaNs` or `Infs` are present, these must be present in corresponding cells
    in both dataframes for the dataframes to evaluate as equal.
    
    Parameters
    ----------
    df1 : :py:class:`pd.DataFrame`
        First dataframe

    df2 : :py:class:`pd.DataFrame`
        Second dataframe
    
    tol : float, optoinal
        Maximum tolerated difference between floating point numbers
        (Default: *1e-8*)
    
    sort_columns : list, optional
        List of column names or indices on which to sort data before comparing values
    
    printer : file-like, optional
        Any logger importing a ``write()`` method
    
    print_verbose : bool, optional
        Print verbose output to stderr (Default: False)
    
    return_verbose : bool, optional
        If *True*, return a list of failure messages 
    
    
    Returns
    -------
    bool
        `True` if dataframes are equal, `False` otherwise
    
    list
        A list of strings explaining how `df1` and `df2` differ. Only
        returned if `return_verbose` is `True`
    """
    df1 = df1.sort(axis=1).sort(axis=0)
    df2 = df2.sort(axis=1).sort(axis=0)

    failures = []
    keys1 = set(df1.columns)
    keys2 = set(df2.columns)
    retval = True
    if keys1 != keys2:
        failstrs = ["Tables contain different columns (unique keys shown below):",
                    "    %s: %s" % (1,", ".join([str(X) for X in keys1 - keys2])),
                    "    %s: %s" % (2,", ".join([str(X) for X in keys2 - keys1]))
                   ] 
        printer.write("\n".join(failstrs))
        failures.extend(failstrs)
        retval = False
    
    if len(df1) != len(df2):
        failstrs = ["Tables contain different numbers of rows:",
                    "    %s: %s" % (1,len(df1)),
                    "    %s: %s" % (2,len(df2))]
        printer.write("\n".join(failstrs))
        failures.extend(failstrs)
        retval = False

    if retval == True:
        if print_verbose == True:    
            printer.write("Sorting data...")
    
        if print_verbose == True:
            printer.write("Testing equality of values by column with numerical tolerance %.3e..." % tol)
            
        unequal = []
        for k in keys1:
            if print_verbose == True:
                printer.write("    Testing column %s..." % k)
            # take out values to get sorted numpy arrays, avoiding weird reordering
            # done behind the scenes by pandas
            s1 = df1[k].values
            s2 = df2[k].values
            if not equal_enough(s1,s2,tol=tol,printer=printer):
                unequal.append(k)
            else:
                if print_verbose == True:
                    printer.write("        %s is identical" % k)
    
        failstr = "Column values are unequal for columns: %s" % (", ".join([str(X) for X in sorted(unequal)]))
        if len(unequal) > 0:
            printer.write(failstr)
            failures.append(failstr)
            retval = False

    if return_verbose == True:
        return retval, failures
    else:
        return retval

def main(argv=sys.argv[1:],verbose=False):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :py:func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    
    verbose : bool, optional
        If `True`, return 
    
    
    Returns
    -------
    int
        `0` if files are identical, `1` otherwise
    
    str
        Only returned if `verbose` is selected. String describing how
        tables are unequal (e.g. which columns failed, et c).
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file1",type=str)
    parser.add_argument("file2",type=str)
    parser.add_argument("-v",dest="verbose",default=False,action="store_true",
                        help="Give verbose output")
    parser.add_argument("--sort_keys",default=None,metavar="key",nargs="+",
                        help="If specified, values will be sorted by the column(s) corresponding to these name or numbers (0-indexed) before comparison")
    parser.add_argument("--exclude",type=str,default=[],nargs="+",metavar="key",
                        help="Key or number (0-indexed) of columns to exclude")
    parser.add_argument("--no_header",default=False,action="store_true",
                        help="If specified, no header row is present. Columns "+\
                             "for all other command-line flags "+\
                             "must be referenced by number (starting at zero) "+\
                             "rather than name, and will be assumed to be in "+\
                             "the same order in both files.")
    parser.add_argument("--tol",type=float,default=1e-8,
                        help="Tolerance by which floats are allowed to differ (Default: 1e-8)")
    
    args = parser.parse_args(argv)

    kwargs = { "sep"       : "\t",
               "index_col" : args.sort_keys,
               "comment"   : "#",
              }
    exclude = args.exclude
    if args.no_header is True:
        if args.sort_keys is not None:
            kwargs["index_col"] = [int(X) for X in args.sort_keys]
        exclude = [int(X) for X in exclude]

        with opener(args.file1) as fh:
            df1 = pd.read_table(fh,header=None,**kwargs)
            
        with opener(args.file2) as fh:
            df2 = pd.read_table(fh,header=None,**kwargs)

    else:
        with opener(args.file1) as fh:
            df1 = pd.read_table(fh,header=0,**kwargs)
            
        with opener(args.file2) as fh:
            df2 = pd.read_table(fh,header=0,**kwargs)
    
    if len(args.exclude) > 0:
        printer.write("Excluding columns %s: " % ", ".join(args.exclude))

    for k in exclude:
        if k in df1:
            df1.pop(k)
        if k in df2:
            df2.pop(k)

    test_result, messages = test_dataframe_equality(df1,df2,printer=printer,print_verbose=True,return_verbose=True,tol=args.tol) #test_3(df1,df2)

    if test_result == True:
        printer.write("Files contain equivalent data.")
        exit_code = 0
    else:
        printer.write("Files non-equivalent.")
        exit_code = 1
    
    if __name__ == "__main__":
        sys.exit(exit_code)
    else:
        if args.verbose == True or verbose == True:
            return exit_code, messages
        else:
            return exit_code
            
if __name__ == "__main__":
    main()
