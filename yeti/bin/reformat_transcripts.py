#!/usr/bin/env python
"""Convert transcripts in BED, BigBed, GTF2, GFF3, or PSL format to BED or GTF2 format"""
import argparse
import inspect
import sys

from yeti.util.scriptlib.argparsers import get_transcripts_from_args,\
                                                      get_annotation_file_parser
from yeti.util.io.openers import opener, argsopener, get_short_name
from yeti.util.io.filters import NameDateWriter
from yeti.util.scriptlib.help_formatters import format_module_docstring

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
    parser.add_argument("--escape",default=False,action="store_true",
                        help="If specified, tokens in column 9 will be URL-escaped (default: False)")
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
                fout.write(transcript.as_gtf(escape=args.escape))
            elif args.output_format == "BED":
                fout.write(transcript.as_bed())
            if c % 1000 == 1:
                printer.write("Processed %s transcripts..." % c)
            c += 1
    
    printer.write("Processed %s transcripts total." % c)
    printer.write("Done.")

if __name__ == "__main__":
    main()
