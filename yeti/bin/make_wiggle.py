#!/usr/bin/env python
"""Converts alignments in bowtie or BAM format to BEDGraph or wiggle files for genome browsing.
Optionally applies read mapping transformations (e.g. for P-site mapping)
and/or normalization.

Output files
------------
Because `Wiggle`_ and `bedGraph`_ files are unstranded, two files are created:

    ${OUTBASE}_fw.wig
        Counts at each position for the plus strand of each chromosome
    
    ${OUTBASE}_rc.wig
        Counts at each position for the minus strand of each chromosome


See Also
--------
:py:mod:`~yeti.genomics.genome_array` : explanations of mapping transformations and why they can be useful
"""
__author__ = "joshua"
__date__ = "2011-03-18"

from yeti.util.scriptlib.argparsers import get_alignment_file_parser,\
                                                      get_genome_array_from_args
from yeti.util.io.filters import NameDateWriter
from yeti.util.io.openers import get_short_name, argsopener
from yeti.util.services.colors import get_rgb255_from_str
from yeti.util.scriptlib.help_formatters import format_module_docstring
import inspect
import sys

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


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
    import argparse
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[get_alignment_file_parser()])
    parser.add_argument("-o","--out",dest="output_file",type=str,required=True,
                        metavar="FILENAME",
                        help="Name stub for output files")
    parser.add_argument("--window_size",default=100000,metavar="N",type=int,
                        help="Size of nucleotides to fetch at once for export. "+\
                             "Large values are faster but require more memory "+\
                             "(default: 100000)")

    track_opts = parser.add_argument_group(title="track options")
    track_opts.add_argument("--color",type=str,default=None,
                        help="An RGB hex string (#NNNNNN, N in [0-9,A-F]) specifying \
                              the track color.")
    track_opts.add_argument("-t","--track_name",dest="track_name",type=str,
                        help="Name to give browser track",
                        default=None)
    track_opts.add_argument("--output_format",choices=("bedgraph","variable_step"),
                        default="bedgraph",
                        help="Format of output file (default: bedgraph)")

    args = parser.parse_args(argv)
    gnd  = get_genome_array_from_args(args,printer=printer)
    
    if args.track_name is None:
        name = args.output_file
    else:
        name = args.track_name
    
    if args.color is not None:
        fw_color = rc_color = str(get_rgb255_from_str(args.color))[1:-1].replace(" ","")
    else:
        fw_color = rc_color = "0,0,0"
    
    if args.output_format == "bedgraph":
        outfn = gnd.to_bedgraph
    elif args.output_format == "variable_step":
        outfn = gnd.to_variable_step

    with argsopener("%s_fw.wig" % args.output_file,args,"w") as fw_out:
        printer.write("Writing forward strand wiggle...")
        outfn(fw_out,"%s_fw" % name,"+",window_size=args.window_size,color=fw_color)
        fw_out.close()

    with argsopener("%s_rc.wig" % args.output_file,args,"w") as rc_out:
        printer.write("Writing reverse strand wiggle...")
        outfn(rc_out,"%s_rc" % name,"-",window_size=args.window_size,color=rc_color)
        rc_out.close()
    
    printer.write("Done!")


if __name__ == "__main__":
    main()
