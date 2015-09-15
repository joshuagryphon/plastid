#!/usr/bin/env python
"""Convert transcripts from `BED`_, `BigBed`_, `GTF2`_, `GFF3`_, or `PSL`_ format
to `BED`_ or `GTF2`_ format.

Notes
-----
GFF3 schemas vary
    Different GFF3s have different schemas of hierarchy. We deal with that here
    by allowing users to supply `transcript_types` and `exon_types`, to indicate
    which sorts of features should be included.

Identity relationships between elements vary between GFF3 files
    Also, different GFF3s specify discontiguous features differently. For example,
    in Flybase, different exons of a transcript will have unique IDs, but will share
    the same "Parent" attribute in column 9 of the GFF. In Wormbase, however, different
    exons of the same transcript will share the same ID. Here, we treat GFFs as if
    they are written in the Flybase style. We may support alternate formats in the future.    

"""
import argparse
import inspect
import warnings
import sys

from plastid.util.scriptlib.argparsers import get_transcripts_from_args,\
                                                      get_annotation_file_parser
from plastid.util.io.openers import argsopener, get_short_name
from plastid.util.io.filters import NameDateWriter
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
    annotation_parser = get_annotation_file_parser()

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     parents=[annotation_parser],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no_escape",default=True,action="store_false",
                        help="If specified and output format is GTF2, special characters in column 9 will be escaped (default: True)")
    parser.add_argument("--output_format",choices=["BED","GTF2"],default="GTF2",
                        help="Format of output file. (default: GTF2)")
    parser.add_argument("outfile",metavar="outfile.gtf",type=str,
                        help="Output file")
    args = parser.parse_args(argv)

    with argsopener(args.outfile,args,"w") as fout:
        c = 0
        transcripts = get_transcripts_from_args(args)
        
        for transcript in transcripts:
            if args.output_format == "GTF2":
                fout.write(transcript.as_gtf(escape=args.no_escape))
            elif args.output_format == "BED":
                fout.write(transcript.as_bed())
            if c % 1000 == 1:
                printer.write("Processed %s transcripts..." % c)
            c += 1
    
    printer.write("Processed %s transcripts total." % c)
    printer.write("Done.")

if __name__ == "__main__":
    main()
