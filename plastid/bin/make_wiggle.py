#!/usr/bin/env python
"""Create :term:`genome browser` tracks from :term:`read alignments`, using
:term:`mapping rules <mapping rule>` to extract the biology of interest
(e.g. ribosomal P-sites, et c) from the alignments.


Output files
------------
Tracks can be output in `wiggle`_ and `bedGraph`_ formats. Because these formats
are unstranded, two files are created:

    OUTBASE_fw.wig
        Counts at each position for the plus/forward strand of each chromosome
    
    OUTBASE_rc.wig
        Counts at each position for the minus/reverse strand of each chromosome

where `OUTBASE` is given by the user.


See also
--------
:doc:`/concepts/mapping_rules`
    Explanations of mapping rules and why they can be useful

:mod:`plastid.genomics.map_factories`
    For lists of mapping rules and their parameters
"""
__author__ = "joshua"
__date__ = "2011-03-18"
import warnings
import inspect
import sys
import argparse

from plastid.util.scriptlib.argparsers import AlignmentParser, BaseParser
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import get_short_name, argsopener
from plastid.plotting.colors import get_rgb255
from plastid.util.scriptlib.help_formatters import format_module_docstring

warnings.simplefilter("once")
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
    ap = AlignmentParser()
    bp = BaseParser()
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[bp.get_parser(),ap.get_parser()])
    parser.add_argument("-o","--out",dest="outbase",type=str,required=True,
                        metavar="FILENAME",
                        help="Base name for output files")
    parser.add_argument("--window_size",default=100000,metavar="N",type=int,
                        help="Size of nucleotides to fetch at once for export. "+\
                             "Large values are faster but require more memory "+\
                             "(Default: 100000)")

    track_opts = parser.add_argument_group(title="Browser track options")
    track_opts.add_argument("--color",type=str,default=None,
                        help="An RGB hex string (`'#NNNNNN'`, `N` in `[0-9,A-F]`) specifying \
                              the track color.")
    track_opts.add_argument("-t","--track_name",dest="track_name",type=str,
                        help="Name to give browser track",
                        default=None)
    track_opts.add_argument("--output_format",choices=("bedgraph","variable_step"),
                        default="bedgraph",
                        help="Format of output file (Default: bedgraph)")

    args = parser.parse_args(argv)
    gnd  = ap.get_genome_array_from_args(args,printer=printer)
    bp.get_base_ops_from_args(args)
    
    if args.track_name is None:
        name = args.outbase
    else:
        name = args.track_name
    
    if args.color is not None:
        fw_color = rc_color = "%s,%s,%s" % tuple(get_rgb255(args.color))
    else:
        fw_color = rc_color = "0,0,0"
    
    if args.output_format == "bedgraph":
        outfn = gnd.to_bedgraph
    elif args.output_format == "variable_step":
        outfn = gnd.to_variable_step

    track_fw = "%s_fw.wig" % args.outbase
    track_rc = "%s_rc.wig" % args.outbase

    with argsopener(track_fw,args,"w") as fw_out:
        printer.write("Writing forward strand track to %s ..." % track_fw)
        outfn(fw_out,"%s_fw" % name,"+",window_size=args.window_size,color=fw_color,
                printer=printer)
        fw_out.close()

    with argsopener(track_rc,args,"w") as rc_out:
        printer.write("Writing reverse strand track to %s ..." % track_rc)
        outfn(rc_out,"%s_rc" % name,"-",window_size=args.window_size,color=rc_color,
                printer=printer)
        rc_out.close()
    
    printer.write("Done!")


if __name__ == "__main__":
    main()
