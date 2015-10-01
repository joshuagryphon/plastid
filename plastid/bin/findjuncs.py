#!/usr/bin/env python 
"""Identify all the unique splice junctions in one or more transcript annotations,
and output these as a `BED`_ file with one splice junction per line. Optionally,
also export junctions as a `Tophat`_ ``.juncs`` file.

If a splice junction is multiply annotated (e.g. used by more than one transcript),
only the first occurrence of the junction will be reported. Scores, if present,
are exported unaltered in `BED`_ output. Examples:

 .. code-block:: shell

    # identify splice junctions from a transcript annotation supplied in GTF2
    # creates output file 'annotation.bed'
    $ findjuncs my_annotation --annotation_format GTF2 \\
                --annotation_files transcripts.gtf
    
    # merge unique annotations from annotation.bed and newly_discovered.bed,
    #
    # export only unique junctions to 'merged_unique.bed'
    $ findjuncs merged_unique --annotation_format BED \\
                --annotation_files annotation.bed newly_discovered.bed


See also
--------
:py:mod:`plastid.bin.slidejuncs`
    Script that makes richer comparisons between discovered and annotated
    junctions, using genomic sequence and :py:mod:`plastid.bin.crossmap`
    results to classify junctions
"""
"""
From the Tophat specification for `.juncs` files:

    Junctions are specified one per line, in a tab-delimited format. Records
    look like::
    
        <chrom> <left> <right> <+/->
    
    left and right are zero-based coordinates, and specify the last character of
    the left sequenced to be spliced to the first character of the right sequence,
    inclusive. That is, the last and the first positions of the flanking exons.
    Users can convert junctions.bed (one of the TopHat outputs) to this format
    using bed_to_juncs < junctions.bed > new_list.juncs where bed_to_juncs can
    be found under the same folder as tophat

See http://ccb.jhu.edu/software/tophat/index.shtml for more information.
"""
import sys
import argparse
import inspect
import warnings
from plastid.genomics.roitools import SegmentChain
from plastid.util.scriptlib.argparsers import get_annotation_file_parser, get_segmentchains_from_args
from plastid.util.io.filters import NameDateWriter
from plastid.util.io.openers import argsopener, get_short_name
from plastid.util.scriptlib.help_formatters import format_module_docstring

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

_ANNOTATION_INPUT_CHOICES = ["BED","BigBed","GTF2","GFF3","PSL"]
_ANNOTATION_DISABLED = ["add_three","annotation_file"]

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
    annotation_file_parser = get_annotation_file_parser(input_choices=_ANNOTATION_INPUT_CHOICES)
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[annotation_file_parser])
    parser.add_argument("--export_tophat",default=False,action="store_true",
                         help="Export tophat `.juncs` file in addition to BED output")
    parser.add_argument("outbase",type=str,help="Basename for output files")

    args = parser.parse_args(argv)
    transcripts = get_segmentchains_from_args(args,printer=printer)
    
    with argsopener("%s.bed" % args.outbase,args,"w") as bed_out:
        if args.export_tophat == True:
            tophat_out = open("%s.juncs" % args.outbase,"w")
    
        printer.write("params: " +" ".join(argv))
        printer.write("Detecting & comparing junctions...")
        ex_pairs = {}
        
        c = 0
        u = 0
        for ivc in transcripts:
            if len(ivc) > 1: # if multi-exon
                chrom = ivc.chrom
                strand = ivc.strand
                for i in range(0,len(ivc)-1):
                    seg1 = ivc[i]
                    seg2 = ivc[i+1]
                    if c % 1000 == 0 and c > 0:
                        printer.write("Processed %s junctions. Found %s unique..." % (c,u) )
                    c+=1
                    key = (chrom,seg1.end,seg2.start,strand)
                    if key not in ex_pairs.keys():
                        u += 1
                        new_ivc = SegmentChain(seg1,seg2)
                        ex_pairs[key] = new_ivc
                        bed_out.write(new_ivc.as_bed())
                        if args.export_tophat == True:
                            my_junc = (chrom,seg1.end-1,seg2.start,strand)
                            tophat_out.write("%s\t%s\t%s\t%s\n" % my_junc)
        del seg1
        del seg2
        del ivc
    
        printer.write("Processed %s total junctions. Found %s unique." % (c,u) )
    
        bed_out.close()
        if args.export_tophat == True:
            tophat_out.close()

    printer.write("Done.")


if __name__ == "__main__":
    main()
