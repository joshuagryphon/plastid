#!/usr/bin/env python
"""Test suite for :py:mod:`plastid.bin.test_table_equality`"""
import os
import shlex
import numpy
import pandas as pd
from random import shuffle
from tempfile import NamedTemporaryFile
from pkg_resources import cleanup_resources
from nose.plugins.attrib import attr
from nose.tools import assert_equal
from plastid.util.io.filters import CommentWriter
from plastid.util.io.openers import NullWriter
from plastid.genomics.roitools import GenomicSegment
from plastid.bin.test_table_equality import main
from plastid.util.services.decorators import catch_stderr



size = 5000

#capture output
run_main = catch_stderr()(main)

@attr(test="functional")
def test_exit_status():

    # define columns
    cols = {
        "intA"   : numpy.random.randint(0,high=2**16,size=size),
        "intB"   : numpy.random.randint(-10,high=20,size=size),
        "idxA"   : numpy.arange(size),
        "chrA"   : numpy.array([chr(65+(X%(91-65))) for X in range(size)]),
        "strA"   : numpy.array([str(GenomicSegment("chrA",X,X+500,"+")) for X in range(size)]),
        "strB"   : numpy.array([str(GenomicSegment("chrB",X/2,X/2+500,"-")) for X in range(size)]),
        "floatA" : 10*numpy.random.randn(size) + 500,
        "floatB" : (10**-5)*numpy.random.random(size),
        "objA"   : numpy.tile(None,5000),
                "objB"   : numpy.array([GenomicSegment("chrC",X,X+Y,"+") for X,Y in zip(range(size),numpy.random.randint(2,high=1000,size=size))]),
          }

    # allocate temp files we will use
    headerfile            = NamedTemporaryFile(delete=False,mode="w")
    headerfile_extra_cols = NamedTemporaryFile(delete=False,mode="w")
    headerfile_extra_cols_diff = NamedTemporaryFile(delete=False,mode="w")
    headerfile_extra_cols_shuffled = NamedTemporaryFile(delete=False,mode="w")
    headerfile_shuffled   = NamedTemporaryFile(delete=False,mode="w")
    headerfile_diff_vals  = NamedTemporaryFile(delete=False,mode="w")

    noheaderfile            = NamedTemporaryFile(delete=False,mode="w")
    noheaderfile_extra_cols = NamedTemporaryFile(delete=False,mode="w")
    noheaderfile_extra_cols_diff = NamedTemporaryFile(delete=False,mode="w")
    noheaderfile_extra_cols_shuffled = NamedTemporaryFile(delete=False,mode="w")
    noheaderfile_shuffled   = NamedTemporaryFile(delete=False,mode="w")
    noheaderfile_diff_vals  = NamedTemporaryFile(delete=False,mode="w")

    # write values
    keyorder = ["idxA"] + sorted(list(set(cols.keys()) - { "idxA" }))

    table1 = pd.DataFrame(cols)
    table1.to_csv(headerfile,index=False,header=True,sep="\t")
    table1.to_csv(noheaderfile,index=False,header=False,sep="\t",
columns=keyorder)
    headerfile.close()
    noheaderfile.close()

    table1["extra"] = 2**7 * numpy.random.random(size=size)
    table1.to_csv(headerfile_extra_cols,index=False,header=True,sep="\t")
    table1.to_csv(noheaderfile_extra_cols,index=False,header=False,sep="\t",

                   columns=["extra"]+keyorder)
    headerfile_extra_cols.close()
    noheaderfile_extra_cols.close()

    table1["extra"] += 10**-4 * numpy.random.random(size=size)
    table1.to_csv(headerfile_extra_cols_diff,index=False,header=True,sep="\t")
    table1.to_csv(noheaderfile_extra_cols_diff,index=False,header=False,sep="\t",

                   columns=["extra"]+keyorder)
    headerfile_extra_cols_diff.close()
    noheaderfile_extra_cols_diff.close()

    shufidx = numpy.arange(size)
    shuffle(shufidx)
    table2 = pd.DataFrame({ K : V[shufidx] for K,V in cols.items()})
    table2.to_csv(headerfile_shuffled,index=False,header=True,sep="\t")
    table2.to_csv(noheaderfile_shuffled,index=False,header=False,sep="\t",
columns=keyorder)
    headerfile_shuffled.close()
    noheaderfile_shuffled.close()

    table2["extra"] = table1["extra"][shufidx]
    table2.to_csv(headerfile_extra_cols_shuffled,index=False,header=True,sep="\t")
    table2.to_csv(noheaderfile_extra_cols_shuffled,
                  index=False,header=False,sep="\t",
                   columns=["extra"]+keyorder)
    headerfile_extra_cols_shuffled.close()
    noheaderfile_extra_cols_shuffled.close()

    # Define tests, as tuples of:
    #   -Test name/description
    #   -Command-line arguments to pass to :py:mod:`plastid.bin.test_table_equality`
    #   -Expected exit code/returns status for :py:func:`main`
    tests = [
        ("same",
            "%s %s" % (headerfile.name,headerfile.name),
            0),
        ("diff_column_names",
            "%s %s" % (headerfile.name,headerfile_extra_cols.name),
            1),
        ("extra_column_names_ignored",
            "%s %s --exclude extra" % (headerfile.name,headerfile_extra_cols.name),
            0),
        ("shuffled_rows",
            "%s %s" % (headerfile.name,headerfile_shuffled.name),
            1),
        ("shuffled_rows_name_sort",
            "%s %s --sort_keys idxA" % (headerfile.name,headerfile_shuffled.name),
            0),
        ("shuffled_rows_multi_name_sort",
            "%s %s --sort_keys strB chrA" % (headerfile.name,headerfile_shuffled.name),
            0),
        ("same_column_names_diff_values",
            "%s %s" % (headerfile_extra_cols.name,headerfile_extra_cols_diff.name),
            1),
        ("same_column_names_diff_values_tol",
            "%s %s --tol 0.01" % (headerfile_extra_cols.name,headerfile_extra_cols_diff.name),
            0),
        ("same_column_names_diff_values_ignored",
            "%s %s --exclude extra" % (headerfile_extra_cols.name,headerfile_extra_cols_diff.name),
            0),
        ("shuffled_rows_extra_columns_ignored",
            "%s %s --exclude extra" % (headerfile.name,headerfile_extra_cols_shuffled.name),
            1),
        ("shuffled_rows_extra_columns_ignored_name_sort",
            "%s %s --exclude extra --sort_keys idxA" % (headerfile.name,headerfile_extra_cols_shuffled.name),
            0),

        ("noheader_same",
            "%s %s --no_header" % (noheaderfile.name,noheaderfile.name),
            0),
        ("noheader_extra_columns",
            "%s %s --no_header" % (noheaderfile.name,noheaderfile_extra_cols.name),
            1),
        ("noheader_shuffled_rows",
            "%s %s --no_header" % (noheaderfile.name,noheaderfile_shuffled.name),
            1),
        ("noheader_shuffled_rows_int_sort",
            "%s %s --no_header --sort_keys 0" % (noheaderfile.name,noheaderfile_shuffled.name),
            0),
        ("noheader_diff_values",
            "%s %s --no_header" % (noheaderfile_extra_cols.name,noheaderfile_extra_cols_diff.name),
            1),
        ("noheader_diff_values_tol",
            "%s %s --no_header --tol 0.01" % (noheaderfile_extra_cols.name,noheaderfile_extra_cols_diff.name),
            0),
        ("no_header_diff_values_ignored",
            "%s %s --no_header --exclude 0" % (noheaderfile_extra_cols.name,noheaderfile_extra_cols_diff.name),
            0),
        ("no_header_shuffled_rows_extra_columns_ignored_int_sort",
            "%s %s --no_header --exclude 0 --sort_keys 1" % (noheaderfile_extra_cols.name,noheaderfile_extra_cols_shuffled.name),
            0)
    ]
    """ Tests to conduct, as tuples of:

        - Test name/description
        - Command-line arguments to pass to :py:mod:`plastid.bin.test_table_equality`
        - Expected exit code/returns status for :py:func:`main`
    """
    for test_name, argstr, expected_exit in tests:
        yield check_exit_status, test_name, argstr, expected_exit

    # clean up
    os.unlink(headerfile.name            )
    os.unlink(headerfile_extra_cols.name )
    os.unlink(headerfile_extra_cols_diff.name )
    os.unlink(headerfile_extra_cols_shuffled.name )
    os.unlink(headerfile_shuffled.name   )
    os.unlink(headerfile_diff_vals.name  )

    os.unlink(noheaderfile.name            )
    os.unlink(noheaderfile_extra_cols.name )
    os.unlink(noheaderfile_extra_cols_diff.name )
    os.unlink(noheaderfile_extra_cols_shuffled.name )
    os.unlink(noheaderfile_shuffled.name   )
    os.unlink(noheaderfile_diff_vals.name  )
    cleanup_resources()

def check_exit_status(test_name,argstr,expected_exit):
    """Verify exit status from :py:mod:`~plastid.bin.test_table_equality`

    Parameters
    ----------
    test_name : str
        Name of test (for reporting purposes)

    argstr : str
        Command-line arguments to pass to :py:mod:`~plastid.bin.test_table_equality`

    expected_exit : int
        Expected exit status
    """
    found_exit, failures = run_main(shlex.split(argstr),verbose=True)
    msg = "    \n".join(failures)
    assert_equal(found_exit,expected_exit,
                 "Failed test %s. Got exit code %s, expected %s.\n    %s" % \
                         (test_name,found_exit,expected_exit,msg))

