#!/usr/bin/env python
"""Compare splice junctions discovered in a dataset to those in an annotation
of known splice junctions, amending misplaced junctions, and identifying 
junctions that fall within repetitive areas of the genome.

Known splice junctions can be misidentified as novel or non-canonical junctions
when intronic sequence immediately downstream of the fiveprime splice site
exactly matches the exonic sequence immediately downstream of the threeprime
splice site. In fact, the junction point could appear anywhere in this
locally-repeated region with equal support from sequencing data. For example,
suppose we have a splice junction as follows::

        
                Exon 1 [0,6)            Intron                                  Exon 2 [16,24)
                ---------------------   --------------------------------------  ------------------------------
    Sequence    G   C   T   C   T   A   C   T   A   G   N   N   N   C   T   A   C   T   A   G   A   T   G   G
    Position    0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23
    Repeated                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    
In this case, the splice junction could be moved 3 bases to the left, or
four bases to the right, without losing consistency with the sequence of
any cDNA or read alignment covering the junction.
   
To identify this and other causes of false positive splice junction calls,
the following operations are performed on each query junction:

    1.  If a :term:`mask file` from |crossmap|
        is provided, junctions in which one or more of the 5' and 3' splice
        sites appear in a repetitive region of the genome are flagged as
        non-informative and written to a separate file. 

    2.  For remaining splice junctions, the extent of locally repeated nucleotide
        sequence, if any, surrounding the query junction's splice donor and
        acceptor sites, are determined in both the 5' and 3' directions.
        
        This is the maximum window (*equal-support region*) across which the
        actual splice junction could be moved without reducing sequence support.
    
    3.  If there is one or more known splice junctions in this region, the
        query junction is assumed to match these, and the known junctions are
        reported rather than the query. 
    
    4.  If (3) is not satisfied, and the query junction is a canonical splice
        junction, it is reported as is.

    5.  If (3) is not satisfied, and the query junction represents a non-canonical
        splice junction, the program determines if one or more canonical splice
        junctions is present in the equal-support region. If so, these canonical
        splice junction are reported rather than the query junction.
        
    6.  If (5) is not satisfied, the non-canonical query junction is reported as-is.

Output files
------------
The following files are written, where `OUTBASE` is a string supplied by the
user. Scores of splice junctions, if present in the input, are ignored.
Each record in each `BED`_ file represents a single exon-exon junction, rather
than a transcript:

    OUTBASE_repetitive.bed
        Splice junctions in which one or more of the splice sites lands
        in a repetitive/degenerate region of the genome, which gives rise to
        mapping ambiguities (step 1 above)
    
    OUTBASE_shifted_known.bed
        The result of shifting query splice junctions to known splice junctions
        with equal sequence support (step 3 above)
    
    OUTBASE_shifted_canonical.bed
        The result of shifting non-canonical query splice junctions to canonical
        splice junctions with equal sequence support (step 5 above)
    
    OUTBASE_untouched.bed
        Query junctions reported without changes (steps 4 and 6 above)

where `OUTBASE` is given by the user.

"""
__date__ =  "Aug 23, 2011"
__author__ = "joshua"

import sys
import argparse
from plastid.util.io.openers import opener, argsopener, get_short_name
from plastid.util.io.filters import CommentReader, NameDateWriter
from plastid.genomics.splicing import get_junction_tuple
from plastid.genomics.roitools import GenomicSegment, SegmentChain
from plastid.genomics.genome_hash import GenomeHash
from plastid.readers.bed import BED_Reader
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.scriptlib.argparsers import (MaskParser, SequenceParser, BaseParser)
import inspect
import warnings
warnings.simplefilter("once")

printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


#===============================================================================
# INDEX : helper functions that classify/manipulate splice junctions 
#===============================================================================

def find_match_range(seg,genome,maxslide):
    """Find maximum distance over which a splice junction can be moved up-
    or down-stream without reducing sequencing support for that junction.
    
    In other words, find locally repeated sequences surrounding exon-intron
    or intron-exon boundaries that can cause splice junction mapping to be
    ambiguous, due to identical and repeated sequence.

    Parameters
    ----------
    seg : |SegmentChain|
         A two-exon fragment representing a query splice junction
         
    genome : dict
        dict mapping chromosome names to :py:class:`Bio.SeqRecord.SeqRecord` s
        
    maxslide : int
        Maximum number of nucleotides from the boundary over which to check
        for extent of repeated sequence
    
    Returns
    -------
    minus_range : int
        Maximum number of nucleotides splice junction point could be moved 
        to the left without reducing sequence support for the junction
        
    plus_range : int
        Maximum number of nucleotides splice junction point could be moved 
        to the right without reducing sequence support for the junction
    """
    iv1,iv2 = seg[0],seg[1]
    chrom = seg.chrom
    
    plus_range  = 0
    minus_range = 0
    check_plus  = True
    check_minus = True
    
    # Discovery of match range. Suppose we have a splice junction as follows:
    #     
    #        Exon 1 [0,6)            Intron                                  Exon 2 [16,24)
    #        ---------------------   --------------------------------------  ------------------------------
    # Seq    G   C   T   C   T   A   C   T   A   G   N   N   N   C   T   A   C   T   A   G   A   T   G   G
    # Pos    0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23
    # 
    # 
    # Exon 1 start: 0
    # Exon 1 end:   6
    # Exon 2 start: 16
    # Exon 2 end:   24
    # 
    # Match range: -3, 4
    # 
    # Then seq[ex1.end    ] == seq[6] == C == seq[16] == seq[ex2.start    ]
    # ---------------------------------------------------------------------
    # Then seq[ex1.end - 3] == seq[3] == C == seq[13] == seq[ex2.start - 3]
    # Then seq[ex1.end - 2] == seq[4] == T == seq[14] == seq[ex2.start - 2]
    # Then seq[ex1.end - 1] == seq[5] == A == seq[15] == seq[ex2.start - 1]
    # Then seq[ex1.end + 0] == seq[6] == C == seq[16] == seq[ex2.start + 0]
    # Then seq[ex1.end + 1] == seq[7] == T == seq[17] == seq[ex2.start + 1]
    # Then seq[ex1.end + 2] == seq[8] == A == seq[18] == seq[ex2.start + 2]
    # Then seq[ex1.end + 3] == seq[9] == G == seq[19] == seq[ex2.start + 3]
    # 
    # We must therefore search for known, repetitive, or canonical splice junctions
    # by iterating our offset ``i`` over [-3,3] = ``range(minus_range, plus_range+1)`` 

    for i in range(maxslide+1):
        if check_plus is True and\
        genome[chrom][iv1.end + i] == genome[chrom][iv2.start + i]:
            plus_range = i+1 #+=1
        else:
            check_plus = False
            
        if check_minus is True and i > 0:
            if genome[chrom][iv1.end - i] == genome[chrom][iv2.start - i]:
                minus_range = -i #-=1
            else:
                check_minus = False
    return minus_range, plus_range            

def find_known_in_range(query_junc,minus_range,plus_range,knownjunctions):
    """Find any known splice junctions within in `minus_range...plus_range`
    of `query_junc`
    
    To be classified as within the range, the boundaries of a known
    junction must be:
    
      1. within `minus_range...plus_range` of the boundaries
         of the the discovered junction.
            
      2. separated by a nucleotide distance equal to the distance 
         separating the junction in `query_junc`.
        
      3. on the same chromosome and strand.
    
    
    Parameters
    ----------
    query_junc : |SegmentChain|
         A two-exon fragment representing a query splice junction
         
    minus_range : int <= 0
        Maximum number of nucleotides splice junction could be moved 
        to the left without reducing sequence support for the junction
        see :py:func:`find_match_range`
        
    plus_range : int >= 0
        Maximum number of nucleotides splice junction could be moved 
        to the right without reducing sequence support for the junction
        see :py:func:`find_match_range`
    
    knownjunctions : list of |SegmentChains|
        known splice junctions
        
    Returns
    -------
    list
        List of |SegmentChains| representing known splice junctions in
        `minus_range...plus_range` of `query_junc`
    """
    exact_match = None
    ltmp = []
    iv1,iv2  = query_junc[0], query_junc[1]
    iv1end   = iv1.end
    iv2start = iv2.start
    qstrand  = query_junc.strand
    qchrom   = query_junc.chrom

    for ivc in knownjunctions:
        if qchrom != ivc.chrom:
            continue
        elif qstrand != ivc.strand:
            continue
        elif query_junc == ivc:
            exact_match = query_junc
            break
        kiv1, kiv2 = ivc[0],ivc[1]
        kend = kiv1.end
        kstart = kiv2.start
        if  kend   in set(range(iv1end+minus_range,  iv1end  +plus_range+1))\
        and kstart in set(range(iv2start+minus_range,iv2start+plus_range+1)):
            if kstart - kend == iv2start - iv1end:
                ltmp.append(ivc)
    return [exact_match] if exact_match is not None else ltmp
        
def find_canonicals_in_range(query_junc,minus_range,plus_range,genome,canonicals):
    """Find any canonical splice junctions within in `minus_range...plus_range`
    of `query_junc`
    
    To be classified as within the range, the boundaries of the canonical
    junction must be:
    
      1. within `minus_range...plus_range` of the boundaries
         of the the discovered junction.
            
      2. separated by a nucleotide distance equal to the distance 
         separating the junction in `query_junc`.
        
      3. On the same chromosome and strand.
        
    
    Parameters
    ----------
    query_junc : |SegmentChain|
         A two-exon fragment representing a query splice junction
         
    minus_range : int <= 0
        Maximum number of nucleotides splice junction could be moved 
        to the left without reducing sequence support for the junction
        see :py:func:`find_match_range`
        
    plus_range : int >= 0
        Maximum number of nucleotides splice junction could be moved 
        to the right without reducing sequence support for the junction
        see :py:func:`find_match_range`
        
    genome : dict
        dict mapping chromosome names to :py:class:`Bio.SeqRecord.SeqRecord` s
        
    canonicals : list
        dinucleotide sequences to consider as canonical splice sites,
        as a list of tuples. e.g. `[("GT","AG"), ("GC","AG")]`

        
    Returns
    -------
    list
        List of |SegmentChains| representing canonical splice junctions in
        `minus_range...plus_range` of `query_junc`
    """    
    ltmp = []
    chrom  = query_junc.chrom
    strand = query_junc.strand
    iv1,iv2 = query_junc[0], query_junc[1]

    iv1start, iv1end = iv1.start, iv1.end
    iv2start, iv2end = iv2.start, iv2.end
    for i in range(minus_range,plus_range+1):
        for pair in canonicals:
            if str(genome[chrom][iv1end + i:iv1end + i + 2].seq) == pair[0]\
            and str(genome[chrom][iv2start - 2 + i:iv2start + i].seq) == pair[1]:
                new_iv1 = GenomicSegment(chrom,
                                          iv1start,
                                          iv1end + i,
                                          strand)
                new_iv2 = GenomicSegment(chrom,
                                          iv2start + i,
                                          iv2end,
                                          strand)
                ltmp.append(SegmentChain(new_iv1,new_iv2))
    
    return ltmp            

def covered_by_repetitive(query_junc,minus_range,plus_range,cross_hash):
    """Determine whether one or both ends of a splice site overlap with
    a repetitive area of the genome.
    
    Parameters
    ----------
    query_junc : |SegmentChain|
         A two-exon fragment representing a query splice junction
    
    minus_range : int <= 0
        Maximum number of nucleotides splice junction could be moved 
        to the left without reducing sequence support for the junction
        see :py:func:`find_match_range`
        
    plus_range : int >= 0
        Maximum number of nucleotides splice junction could be moved 
        to the right without reducing sequence support for the junction
        see :py:func:`find_match_range`
    
    cross_hash : |GenomeHash|
        |GenomeHash| of 1-length features denoting repetitive regions of the genome
        
    
    Returns
    -------
    bool
        `True` if any of the genomic positions within `minus_range...plus_range`
        of the 5' or 3' splice sites of `query_junc` overlap a repetitive
        region of the genome as annotated by ``cross_hash``.
        Otherwise, `False`
    """
    chrom = query_junc.spanning_segment.chrom
    strand = query_junc.spanning_segment.strand
    qend = query_junc[0].end
    qstart = query_junc[1].start
    fiveprime_splice_area = GenomicSegment(chrom,
                                           qend + minus_range,
                                           qend + plus_range + 1,
                                           strand)
    threeprime_splice_area = GenomicSegment(chrom,
                                            qstart + minus_range,
                                            qstart + plus_range + 1,
                                            strand)
    support_region = SegmentChain(fiveprime_splice_area,threeprime_splice_area)
    return len(cross_hash.get_overlapping_features(support_region)) > 0

#===============================================================================
# INDEX : program body 
#===============================================================================

def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :py:func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """
    sp = SequenceParser()
    mp = MaskParser()
    bp = BaseParser()
    
    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     parents=[bp.get_parser(),sp.get_parser(),mp.get_parser()],
                                     )
    parser.add_argument("--maxslide",type=int,default=10,
                        help="Maximum number of nt to search 5\' and 3\' of intron"+
                             " boundaries (Default: 10)")
    parser.add_argument("--ref",type=str,metavar="ref.bed",default=None,
                        help="Reference file describing known splice junctions")
    parser.add_argument("--slide_canonical",action="store_true",default=False,
                        help="Slide junctions to canonical junctions if present within equal support region")
    parser.add_argument("infile",type=str,metavar="input.bed",
                        help="BED file describing discovered junctions")
    parser.add_argument("outbase",type=str,
                        help="Basename for output files")
    args = parser.parse_args(argv)
    bp.get_base_ops_from_args(args)
    
    printer.write("Opening genome from %s..." % args.sequence_file)
    genome = sp.get_seqdict_from_args(args)
    
    # load crossmap    
    cross_hash = mp.get_genome_hash_from_args(args)

    # load ref junctions
    if args.ref is not None:
        printer.write("Loading reference junctions from %s" % args.ref)
        known_hash = GenomeHash(list(BED_Reader(open(args.ref))),do_copy=False)
    else:
        known_hash = GenomeHash()

    # set up variables    
    canonicals_plus = [("GT","AG"),
                       ("GC","AG")
                      ]
    
    canonicals_minus = [("CT","AC"),
                        ("CT","GC")
                       ]
    
    known_in_range     = 0
    canonical_in_range = 0
    repetitive         = 0
    untouched          = 0
    c = 0
    
    seen_already = []

    outfiles = {
                 "repetitive" : "%s_repetitive.bed" % args.outbase,
                 "known"      : "%s_shifted_known.bed" % args.outbase,
                 "canonical"  : "%s_shifted_canonical.bed" % args.outbase,
                 "untouched"  : "%s_untouched.bed" % args.outbase,
                }
    outfiles = { K : argsopener(V,args,"w") for K,V in outfiles.items() }

    # process data
    printer.write("Opening junctions from %s..." % args.infile)
    for ivc in BED_Reader(CommentReader(opener(args.infile))):
        processed = False
        tup = None

        if c % 1000 == 0 and c > 0:
            printer.write("Processed: %s\tknown: %s\tshifted to canonical: %s\trepetitive: %s\tuntouched: %s" % \
                    (c, known_in_range, canonical_in_range, repetitive, untouched))
                   
        assert len(ivc) == 2
        strand = ivc.strand
        
        minus_range, plus_range = find_match_range(ivc,genome,args.maxslide)
        
        # see if either end of splice junction +- match_range lands in repetitive areas of genome
        if covered_by_repetitive(ivc,minus_range,plus_range,cross_hash):
            repetitive += 1
            outfiles["repetitive"].write(ivc.as_bed())
            processed = True

        # see if one or more known junctions in range
        if processed == False and args.ref is not None:
            # find_known_in_range(query_ivc,minus_range,plus_range,knownjunctions)
            known_juncs = find_known_in_range(ivc,minus_range,plus_range,known_hash.get_nearby_features(ivc))
            if len(known_juncs) > 0:
                known_in_range += 1
                for my_known in known_juncs:
                    tup = get_junction_tuple(my_known)
                    if tup not in seen_already:
                        outfiles["known"].write(my_known.as_bed())
                        seen_already.append(tup)
                    
                processed = True
            
        # see if one or more canonical junctions in range
        if processed == False and args.slide_canonical == True:
            canonicals = canonicals_plus if strand == "+" else canonicals_minus
            #find_canonicals_in_range(query_ivc,minus_range,plus_range,genome,canonicals)
            canonical_juncs = find_canonicals_in_range(ivc,minus_range,plus_range,genome,canonicals)
            if len(canonical_juncs) > 0:
                canonical_in_range += 1
                for can in canonical_juncs:
                    tup = get_junction_tuple(can)
                    if tup not in seen_already:
                        outfiles["canonical"].write(can.as_bed())
                        seen_already.append(tup)

                processed = True
                    
        if processed == False:
            outfiles["untouched"].write(ivc.as_bed())
            untouched += 1
            
        c += 1

    # save output
    printer.write("Totals: %s\tknown: %s\tshifted to canonical: %s\trepetitive: %s\tuntouched: %s" % \
            (c, known_in_range, canonical_in_range, repetitive, untouched))    

    for v in outfiles.values():
        v.close()
    
    printer.write("Done.")


if __name__ == "__main__":
    main()
