#!/usr/bin/env python
"""Retrieve masked vectors of mapped sequencing reads for regions of interest (ROIs).
Vectors are saved as individual files, one for each region of interest,
in a user-specified output folder. If a mask annotation -- e.g. from  
:py:mod:`~yeti.bin.crossmap` -- is provided, masked positions will be
saved as ``nan`` s in output.


Output files
------------
An output folder is created, and the count vector for each ROI is saved as
an individual line-delimited file in the output folder. Each output file is
named after the ROI for which it is created.
"""
import argparse
import inspect
import os
import numpy
import sys

from yeti.util.scriptlib.argparsers import get_alignment_file_parser,\
                                                      get_genome_array_from_args,\
                                                      get_segmentchains_from_args,\
                                                      get_segmentchain_file_parser,\
                                                      get_annotation_file_parser,\
                                                      get_mask_file_parser,\
                                                      get_genome_hash_from_mask_args
from yeti.util.io.openers import opener, get_short_name, argsopener
from yeti.util.io.filters import NameDateWriter
from yeti.util.scriptlib.help_formatters import format_module_docstring


printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

def main(args=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :py:func:`main` is called directly.

        Default: sys.argv[1:] (actual command-line arguments)
    """
    alignment_file_parser  = get_alignment_file_parser()
    annotation_file_parser = get_segmentchain_file_parser()
    mask_file_parser = get_mask_file_parser()
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[alignment_file_parser,
                                              annotation_file_parser,
                                              mask_file_parser])
    
    parser.add_argument("out_folder",type=str,help="Folder in which to put output")
    parser.add_argument("--out_prefix",default="",type=str,
                        help="Prefix to prepend to output files (default: none)")
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
         
        count_vec = tx.get_valid_counts(ga)
        numpy.savetxt(full_filename,count_vec,fmt=args.format)

if __name__ == "__main__":
    main()
