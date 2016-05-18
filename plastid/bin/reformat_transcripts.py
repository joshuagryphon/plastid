#!/usr/bin/env python
"""Convert transcripts from `BED`_, `BigBed`_, `GTF2`_, `GFF3`_, or `PSL`_ format
to `BED`_, :term:`extended BED`, or `GTF2`_ format.

 .. note::

    GFF3 schemas vary
        Different GFF3s have different schemas of hierarchy. By default, we
        assume the ontology used by the Sequence Ontology consortium. Users
        that require a different schema may supply `transcript_types` and
        `exon_types`, to indicate which sorts of features should be included.

    Identity relationships between elements vary between GFF3 files
        GFF3 files can represent discontiguous features using two strategies. In 
        one strategy, the exons of a transcript have unique IDs, but will share
        contain the same parent ID in their same `Parent` attribute in column 9
        of the GFF. In another strategy different exons of the same transcript
        simply share the same ID, and don't define a `Parent`. Here, both schemes
        are accepted, although what happens if they conflict within a single
        transcript is undefined.

"""
import argparse
import inspect
import warnings
import sys
import os
import copy

from plastid.util.scriptlib.argparsers import (AnnotationParser,BaseParser)
from plastid.util.services.exceptions import ArgumentWarning, warn
from plastid.util.io.openers import argsopener, get_short_name
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

MAKE_BIGBED_MESSAGE =\
"""The BED-formatted file `%(outbase)s` has been exported with extended BED
columns. Plastid and other software can detect these automatically if you convert
`%(outbase)s` to a BigBed file using Jim Kent's bedToBigBed utility, as follows:

    $ sort -k1,1 -k2,2n -t"\t" %(outbase)s >%(outbase)s_sorted.bed
    $ bedToBigBed -tab -type=bed12+%(numcols)s -as=fields.as %(outbase)s.bed -extraIndex=name chrom.sizes %(outbase)s.bb

where `chrom.sizes` is a two-column, tab-delimited table of chromosome name & size,
and `fields.as` is an autoSql declaration describing the types & names of the extra
columns. For `%(outbase)s.bed`, `fields.as` could be:

%(autosql)s

Simply edit the type of each custom field to match your data type. Field type
options are:

    int, uint, short, ushort, byte, ubyte, float, char, string, lstring, enum, set

For details, See the autoSql Grammar specification at:

https://github.com/ENCODE-DCC/kentUtils/blob/36d6274459f644d5400843b8fa097b380b8f7867/src/hg/autoSql/autoSql.doc
"""

BED12_RESERVED_NAMES = ["chrom","chromStart","chromEnd","name","score","strand",
                        "thickStart","thickEnd","reserved","blockCount",
                        "blockSizes","chromStarts"]

DEFAULT_AUTOSQL_STR =\
"""
table bigbed_columns "%s columns"
(
    string            chrom;          "Chromosome"
    uint              chromStart;     "chr start"
    uint              chromEnd;       "chr end"
    string            name;           "item name"
    uint              score;          "score"
    char[1]           strand;         "strand"
    uint              thickStart;     "thickstart"
    uint              thickEnd;       "thickend"
    uint              reserved;       "normally itemRgb"
    int               blockCount;     "block count"
    int[blockCount]   blockSizes;     "block sizes"
    int[blockCount]   chromStarts;    "block starts"
%s
)
"""

AUTOSQL_ROW_FMT_STR = '    string            %s;%s"description of custom field contents"' 


def fix_name(inp,names_used):
    """Append a number if an autoSql field name is duplicated.
    """
    name = inp
    i = 2
    while name in names_used:
        name = "%s%s" % (inp,i)
        i += 1
    
    names_used.append(name)
    return name
    
#: TODO: no functional test for --extra_columns
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
    ap = AnnotationParser()
    bp = BaseParser()
    annotation_parser = ap.get_parser()
    base_parser = bp.get_parser()

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     parents=[base_parser,annotation_parser],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--no_escape",default=True,action="store_false",
                        help="If specified and output format is GTF2, special characters in column 9 will be escaped (default: True)")
    parser.add_argument("--output_format",choices=["BED","GTF2"],default="GTF2",
                        help="Format of output file. (default: GTF2)")
    parser.add_argument("--extra_columns",nargs="+",default=[],type=str,
                        help="Attributes (e.g. 'gene_id' to output as extra columns in extended BED format (BED output only).")
    parser.add_argument("--empty_value",default="na",type=str,
                        help="Value to use of an attribute in `extra_columns` is not defined for a particular record (Default: 'na'")
    parser.add_argument("outfile",metavar="outfile.[ bed | gtf ]",type=str,
                        help="Output file")
    args = parser.parse_args(argv)
    bp.get_base_ops_from_args(args)

    end_message = ""    
    extra_cols = args.extra_columns
    if extra_cols is not None:
        if args.output_format == "BED":
            
            # avoid name clashes
            names_used = copy.copy(BED12_RESERVED_NAMES)
            asql_names = [fix_name(X,names_used) for X in extra_cols]
            autosql_str = "\n".join(AUTOSQL_ROW_FMT_STR % (X," "*max(15-len(X),2)) for X in asql_names)
            
            file_info = {
                "outbase" : args.outfile.replace(".bed","").replace(".gtf",""),
                "numcols" : len(extra_cols),
                "autosql" : DEFAULT_AUTOSQL_STR % (os.path.basename(args.outfile[:-4]),autosql_str),
                         
            }
            end_message = MAKE_BIGBED_MESSAGE % file_info
        else:
            warn("`--extra_columns` is ignored for %s-formatted output." % (args.output_format),ArgumentWarning)
            
            
    with argsopener(args.outfile,args,"w") as fout:
        c = 0
        transcripts = ap.get_transcripts_from_args(args,printer=printer)
        
        for transcript in transcripts:
            if args.output_format == "GTF2":
                fout.write(transcript.as_gtf(escape=args.no_escape))
            elif args.output_format == "BED":
                fout.write(transcript.as_bed(extra_columns=extra_cols,empty_value=args.empty_value))
            if c % 1000 == 1:
                printer.write("Processed %s transcripts ..." % c)
            c += 1
    
    printer.write("Processed %s transcripts total." % c)
    printer.write("Done.")
    print(end_message)
            
if __name__ == "__main__":
    main()
