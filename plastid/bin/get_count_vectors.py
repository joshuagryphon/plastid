#!/usr/bin/env python
"""Fetch vectors of :term:`counts` at each nucleotide position in one or more
regions of interest (ROIs).


Output files
------------
Vectors are saved as individual line-delimited files -- one position per line --
in a user-specified output folder. Each file is named for the ROI to which it
corresponds. If a :term:`mask file` -- e.g. from  :py:mod:`~plastid.bin.crossmap`
-- is provided, masked positions will be have value `nan` in output.
"""
import argparse
import inspect
import os
import warnings
import sys
import numpy

from plastid.util.scriptlib.argparsers import get_alignment_file_parser,\
                                                      get_genome_array_from_args,\
                                                      get_segmentchains_from_args,\
                                                      get_segmentchain_file_parser,\
                                                      get_mask_file_parser,\
                                                      get_genome_hash_from_mask_args
from plastid.util.io.openers import get_short_name
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

def main(args=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """
    alignment_file_parser  = get_alignment_file_parser()
    annotation_file_parser = get_segmentchain_file_parser()
    mask_file_parser = get_mask_file_parser()
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[alignment_file_parser,
                                              annotation_file_parser,
                                              mask_file_parser])
    
    parser.add_argument("out_folder",type=str,help="Folder in which to save output vectors")
    parser.add_argument("--out_prefix",default="",type=str,
                        help="Prefix to prepend to output files (default: no prefix)")
    parser.add_argument("--format",default="%.8f",type=str,
                        help=r"printf-style format string for output (default: '%%.8f')")
    args = parser.parse_args(args)

    # if output folder doesn't exist, create it
    if not os.path.isdir(args.out_folder):
        os.mkdir(args.out_folder)
 
    # parse args
    ga = get_genome_array_from_args(args,printer=printer)
    transcripts = get_segmentchains_from_args(args,printer=printer)
    mask_hash = get_genome_hash_from_mask_args(args,printer=printer)
    
    # evaluate
    for n,tx in enumerate(transcripts):
        if n % 1000 == 0:
            printer.write("Processed %s regions of interest" % n)
        filename      = "%s%s.txt" % (args.out_prefix,tx.get_name())
        full_filename = os.path.join(args.out_folder,filename)
         
        # mask out overlapping masked regions
        overlapping = mask_hash.get_overlapping_features(tx)
        for feature in overlapping:
            tx.add_masks(*feature._segments)
         
        count_vec = tx.get_masked_counts(ga)
        numpy.savetxt(full_filename,count_vec,fmt=args.format)

if __name__ == "__main__":
    main()
