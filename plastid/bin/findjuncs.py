#!/usr/bin/env python 
"""This script identify all the unique splice junctions in one or more transcript
annotations, and exports these as a `BED`_ file with one splice junction per line
Optionally, this script can also export junctions as a `Tophat`_ ``.juncs`` file.

If a splice junction appears multiple times (e.g. used by more than one transcript),
only the first occurrence of the junction will be reported. Scores, if present,
are exported unaltered in `BED`_ output.

Examples:


.. code-block:: shell

   # identify splice junctions from a transcript annotation supplied in GTF2
   # creates output file 'annotation.bed'
   $ findjuncs my_annotation --annotation_format GTF2 \\
               --annotation_files transcripts.gtf
    
   # merge unique annotations from annotation.bed and newly_discovered.bed,
   # export only unique junctions to 'merged_unique.bed'
   $ findjuncs merged_unique --annotation_format BED \\
               --annotation_files annotation.bed newly_discovered.bed


See also
--------
:mod:`plastid.bin.slidejuncs`
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
from plastid.util.scriptlib.argparsers import AnnotationParser, BaseParser
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
    ap = AnnotationParser(input_choices=_ANNOTATION_INPUT_CHOICES)
    annotation_file_parser = ap.get_parser()
    
    bp = BaseParser()
    base_parser = bp.get_parser()
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[base_parser,annotation_file_parser])
    parser.add_argument("--export_tophat",default=False,action="store_true",
                         help="Export tophat `.juncs` file in addition to BED output")
    parser.add_argument("outbase",type=str,help="Basename for output files")

    args = parser.parse_args(argv)
    bp.get_base_ops_from_args(args)
    
    transcripts = ap.get_transcripts_from_args(args,printer=printer,return_type=SegmentChain)
    
    with argsopener("%s.bed" % args.outbase,args,"w") as bed_out:
        if args.export_tophat == True:
            tophat_out = open("%s.juncs" % args.outbase,"w")
    
        printer.write("params: " +" ".join(argv))
        printer.write("Detecting & comparing junctions...")
        ex_pairs = {}
        
        c = 0
        u = 0
        for chain in transcripts:
            if len(chain) > 1: # if multi-exon
                chrom = chain.chrom
                strand = chain.strand
                try:
                    ep = ex_pairs[(chrom,strand)]
                except KeyError:
                    ex_pairs[(chrom,strand)] = []
                    ep = ex_pairs[(chrom,strand)]

                for i in range(0,len(chain)-1):
                    
                    seg1 = chain[i]
                    seg2 = chain[i+1]
                    if c % 1000 == 0 and c > 0:
                        printer.write("Processed %s junctions. Found %s unique..." % (c,u) )
                    c+=1
                    key = (seg1.end,seg2.start)
                        
                    if key not in ep:
                        ep.append(key)
                        u += 1
                        new_chain = SegmentChain(seg1,seg2)
                        bed_out.write(new_chain.as_bed())
                        if args.export_tophat == True:
                            my_junc = (chrom,seg1.end-1,seg2.start,strand)
                            tophat_out.write("%s\t%s\t%s\t%s\n" % my_junc)
                        
                        del new_chain
                    
                    del seg1
                    del seg2
                    
            del chain
    
        printer.write("Processed %s total junctions. Found %s unique." % (c,u) )
    
        bed_out.close()
        if args.export_tophat == True:
            tophat_out.close()

    printer.write("Done.")


if __name__ == "__main__":
    main()
