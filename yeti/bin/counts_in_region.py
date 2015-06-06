#!/usr/bin/env python
"""Computes read counts and densities (in reads per nucleotide and in RPKM)
for arbitrary features in an annotation file.
"""
import argparse
import warnings
import inspect
import sys
import itertools
import numpy

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from yeti.util.io.filters import CommentReader, NameDateWriter
    from yeti.util.io.openers import opener, argsopener, get_short_name, guess_opener
    from yeti.util.scriptlib.argparsers import get_genome_array_from_args,\
                                                      get_segmentchains_from_args,\
                                                      get_genome_hash_from_mask_args,\
                                                      get_alignment_file_parser,\
                                                      get_segmentchain_file_parser,\
                                                      get_mask_file_parser
    from yeti.util.scriptlib.help_formatters import format_module_docstring

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

_DISABLED=["normalize"]
def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :py:func:`main` is called directly.

        Default: sys.argv[1:] (actually command-line arguments)
    """
    annotation_file_parser = get_segmentchain_file_parser()
    alignment_file_parser  = get_alignment_file_parser(disabled=_DISABLED)
    mask_file_parser       = get_mask_file_parser()
    
    parser = argparse.ArgumentParser(format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[alignment_file_parser,
                                              annotation_file_parser,
                                              mask_file_parser],
                                     )
    parser.add_argument("outfile",type=str,help="Output filename")
    args = parser.parse_args(argv)
    gnd = get_genome_array_from_args(args,printer=printer,disabled=_DISABLED)
    
    transcripts = get_segmentchains_from_args(args,printer=printer)
    crossmap = get_genome_hash_from_mask_args(args,printer=printer)
    
    gnd_sum = gnd.sum()
    with argsopener(args.outfile,args,"w") as fout:
        fout.write("## total_dataset_counts: %s\n" % gnd_sum)
        fout.write("#region_name\tregion\tcounts\tcounts_per_nucleotide\trpkm\tlength\n")
        for n,ivc in enumerate(transcripts):
            name = ivc.get_name()
            masks = crossmap.get_overlapping_features(ivc)
            ivc.add_masks(*itertools.chain.from_iterable((X for X in masks)))
            if n % 1000 == 0:
                printer.write("Processed %s regions..." % n)
                
            counts = numpy.nansum(ivc.get_valid_counts(gnd))
            length = ivc.get_valid_length()
            rpnt = numpy.nan if length == 0 else float(counts)/length
            rpkm = numpy.nan if length == 0 else float(counts)/length/gnd_sum * 1000 * 10**6 
            ltmp = [name,
                    str(ivc),
                    "%.8e" % counts,
                    "%.8e" % rpnt,
                    "%.8e" % rpkm,
                    "%d" % length]
            fout.write("%s\n" % "\t".join(ltmp))
    
        fout.close()
        
    printer.write("Processed %s regions total." % n)

    printer.write("Done.")

if __name__ == '__main__':
    main()
