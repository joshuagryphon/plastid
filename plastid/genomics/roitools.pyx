#!/usr/bin/env python
"""This module contains classes for representing and manipulating genomic features.

.. contents::
   :local:
   
   
Summary
-------

Genomic features are represented as |SegmentChains|, which can contain zero or
more continuous spans of the genome (|GenomicSegments|), as well as rich
annotation data. For the specific case of RNA transcripts, a subclass of
|SegmentChain|, called |Transcript| is provided.

Module contents
...............

.. autosummary::
   
   GenomicSegment
   SegmentChain
   Transcript
   positions_to_segments
   add_three_for_stop_codon


Examples
--------

|SegmentChains| may be read directly from annotation files using the readers
in :mod:`plastid.readers`::

    >>> from plastid import *
    >>> chains = list(BED_Reader(open("some_file.bed")))

or constructed from |GenomicSegments|:: 

    >>> seg1 = GenomicSegment("chrA",5,200,"-")
    >>> seg2 = GenomicSegment("chrA",250,300,"-")
    >>> my_chain = SegmentChain(seg1,seg2,ID="some_chain", ... , some_attribute="some_value")

|SegmentChains| contain convenience methods for a number of comman tasks, for
example:

  - converting coordinates between the spliced space of the chain, and the genome::
  
        >>> # get coordinate of 50th position from 5' end
        >>> my_chain.get_genomic_coordinate(50)
        ('chrA', 199, '-')
                
        # get coordinate of 49th position. splicing is taken care of!
        >>> my_chain.get_genomic_coordinate(49)
        ('chrA', 250, '-')
        
        # get coordinate in chain corresponding to genomic coordinate 118
        >>> my_chain.get_segmentchain_coordinate("chrA",118,"-")
        131
        
        # get a subchain containing positions 45-70
        >>> subchain = my_chain.get_subchain(45,70)
        >>> subchain
        <SegmentChain segments=2 bounds=chrA:180-255(-) name=some_chain_subchain>

        # the subchain preserves the discontinuity found in `my_chain`
        >>> subchain.segments
        [<GenomicSegment chrA:180-200 strand='-'>,
         <GenomicSegment chrA:250-255 strand='-'>]

  - fetching :class:`numpy arrays <numpy.ndarray>` of data at each position
    in the chain. The data is assumed to be kept in a
    :mod:`GenomeArray <plastid.genomics.genome_array>`::
    
        >>> ga = BAMGenomeArray(["some_file.bam"],mapping=ThreePrimeMapFactory(offset=15))
        >>> my_chain.get_counts(ga)
        array([843, 854, 153,  86, 462, 359, 290,  38,  38, 758, 342, 299, 430,
               628, 324, 437, 231, 417, 536, 673, 243, 981, 661, 415, 207, 446,
               197, 520, 653, 468, 863,   3, 272, 754, 352, 960, 966, 913, 367,
               ...
               ])
               
  - similarly, fetching spliced sequence, reverse-complemented if necessary
    for minus-strand features. As input, the |SegmentChain| expects a
    dictionary-like object mapping chromosome names to string-like sequences
    (e.g. as in `BioPython`_ or `twobitreader`_)::
    
        >>> seqdict = { "chrA" : "TCTACATA ..." } # some string of chrA sequence
        >>> my_chain.get_sequence(seqdict)
        "ACTGTGTACTGTACGATCGATCGTACGTACGATCGATCGTACGTAGCTAGTCAGCTAGCTAGCTAGCTGA..." 
        
  - testing for overlap, containment, equality with other |SegmentChains|::
  
        >>> other_chain = SegmentChain(GenomicSegment("chrA",200,300,"-"),
        >>>                            GenomicSegment("chrA",800,900,"-"))
                                       
        >>>  my_chain.overlaps(other_chain)
        True
        
        >>> other_chain in my_chain
        False
        
        >>> my_chain in my_chain
        True
        
        >>> my_chain.covers(other_chain)
        False
        
        >>> my_chain == other_chain
        False
        
        >>> my_chain == my_chain
        True

  - export to `BED`_, `GTF2`_, or `GFF3`_::
  
        >>> my_chain.as_bed()
        chrA    5    300    some_chain    0    -    5    5    0,0,0    2    195,50,    0,245,
        
        >>> my_chain.as_gtf()
        chrA    .    exon    6    200    .    -    .    gene_id "gene_some_chain"; transcript_id "some_chain"; some_attribute "some_value"; ID "some_chain";
        chrA    .    exon    251  300    .    -    .    gene_id "gene_some_chain"; transcript_id "some_chain"; some_attribute "some_value"; ID "some_chain";
"""



import re
import copy
import array
import numpy
 
cimport numpy
 
from numpy.ma import MaskedArray as MaskedArray
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
 
from plastid.util.services.exceptions import DataWarning, warn
from plastid.util.services.decorators import deprecated
from plastid.plotting.colors import get_str_from_rgb255, get_str_from_rgb, get_rgb255
from plastid.readers.gff_tokens import make_GFF3_tokens, \
                                       make_GTF2_tokens
 
from plastid.genomics.c_common cimport ExBool, true, false, bool_exception, \
                                       Strand, forward_strand, reverse_strand,\
                                       unstranded, undef_strand

#===============================================================================
# constants
#===============================================================================

cdef hash_template = array.array('l',[]) # c signed long / python int
cdef mask_template = array.array('i',[]) # TODO: change to c unsigned char / python int


# to/from str
igvpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)")
segpat = re.compile(r"([^:]*):([0-9]+)-([0-9]+)\(([+-.])\)")
ivcpat = re.compile(r"([^:]*):([^(]+)\(([+-.])\)")

# compile time constants - __richcmp__ test
DEF LT  = 0
DEF LEQ = 1
DEF EQ  = 2
DEF NEQ = 3
DEF GT  = 4
DEF GEQ = 5

# numeric types
INT    = numpy.int
FLOAT  = numpy.float
DOUBLE = numpy.double
LONG   = numpy.long

ctypedef numpy.int_t    INT_t
ctypedef numpy.float_t  FLOAT_t
ctypedef numpy.double_t DOUBLE_t
ctypedef numpy.long_t   LONG_t

# cdef class SegmentChain
# cdef class Transcript(SegmentChain)

#==============================================================================
# Exported object
#==============================================================================

NullSegment = GenomicSegment("NullChromosome",0,0,"\x00")
"""Placeholder |GenomicSegment| for undefined objects"""

#===============================================================================
# io functions
#===============================================================================

def _get_attr(dict dtmp):
    ltmp=[]
    for k,v in sorted(dtmp.items()):
        ltmp.append("%s='%s'" % (k,v))
    
    return ",".join(ltmp)

def _format_segmentchain(SegmentChain segchain):
    """Formats a |SegmentChain| as a string, which when ``eval`` ed, should reconstruct the |SegmentChain|.
    Used for creating test datasets.
    
    Parameters
    ----------
    segchain : |SegmentChain|
    
    Returns
    -------
    str
    """
    iv_ltmp = []
    for iv in segchain:
        iv_ltmp.append("GenomicSegment('%s',%s,%s,'%s')" % (iv.chrom,iv.start,iv.end,iv.strand))
    
    stmp = "SegmentChain(%s,%s)\n" % (",".join(iv_ltmp),
                                    _get_attr(segchain.attr)
                                                          )
    return stmp

def _format_transcript(Transcript tx):
    """Formats a |Transcript| as a string, which when ``eval`` ed, should reconstruct the |Transcript|.
    Used for creating test datasets.
    
    Parameters
    ----------
    tx : |Transcript|
    
    Returns
    -------
    str
    """
    iv_ltmp = []
    for iv in tx:
        iv_ltmp.append("GenomicSegment('%s',%s,%s,'%s')" % (iv.chrom,iv.start,iv.end,iv.strand))
    
    stmp = "Transcript(%s,ID='%s',cds_genome_start=%s,cds_genome_end=%s)\n" % (",".join(iv_ltmp),
                                                          tx.get_name(),
                                                          tx.cds_genome_start,
                                                          tx.cds_genome_end
                                                          )
    return stmp


#==============================================================================
# Exported helper functions e.g. for sorting or building |GenomicSegments|
#==============================================================================

cpdef list merge_segments(list segments):
    """Merge all overlapping |GenomicSegments| in `segments`, so that all segments
    returned are guaranteed to be sorted and non-overlapping.

    .. note::

        All segments are assumed to be on the same strand and chromosome.


    Parameters
    ----------
    segments : list
        List of |GenomicSegments|, all on the same strand and chromosome

    Returns
    -------
    list
        List of sorted, non-overlapping |GenomicSegments|
    """
    cdef:
        GenomicSegment left
        GenomicSegment right
        str chrom
        str strand
        list new_segments = []
        list ssegments
        int i = 0
        int length = len(segments)

    if length < 2:
        return segments

    ssegments = sorted(segments)
    left = ssegments[0]
    chrom = left.chrom
    strand = left.strand

    while i < length - 1:
        right = ssegments[i+1]
        if right.start > left.end:
            new_segments.append(left)
            left = right
        else:
            left = GenomicSegment(chrom,left.start,max(left.end,right.end),strand)
            
        i += 1

    new_segments.append(left)

    return new_segments

cpdef list positions_to_segments(str chrom, str strand, object positions):
    """Construct |GenomicSegments| from a chromosome name, a strand, and a list of chromosomal positions

    Parameters
    ----------
    chrom : str
        Chromosome name
       
    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
       
    positions : list of integers
        **End-inclusive** list, tuple, or set of positions to include
        in final |GenomicSegment|
           
    Returns
    -------
    list
        List of |GenomicSegments| covering `positions`
    """
    if isinstance(positions,set):
        positions = sorted(positions)
    else:
        positions = sorted(set(positions))
    return positionlist_to_segments(chrom,strand,positions)

cpdef list positionlist_to_segments(str chrom, str strand, list positions):
    """Construct |GenomicSegments| from a chromosome name, a strand, and a list of chromosomal positions

    Parameters
    ----------
    chrom : str
        Chromosome name
       
    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
       
    positions : list of unique integers
        Sorted, **end-inclusive** list of positions to include
        in final |GenomicSegment|
           
    Returns
    -------
    list
        List of |GenomicSegments| covering `positions`


     .. warning::
        
        This function is meant to quickly without excessive type conversions.
        So, the elements `positions` must be **UNIQUE** and **SORTED**. If
        they are not, use :func:`positions_to_segments` instead.
    """
    cdef:
        list segments = []
        long start, last_pos, i

    if len(positions) > 0:
        start = positions[0]
        last_pos = positions[0]
        for i in positions:
            if i == start:
                continue
            if i - last_pos == 1:
                last_pos = i
            else:
                segments.append(GenomicSegment(chrom,start,last_pos+1,strand))
                start = i
                last_pos = i
        segments.append(GenomicSegment(chrom,start,last_pos+1,strand))

    return segments


cpdef Transcript add_three_for_stop_codon(Transcript tx):
    """Extend an annotated CDS region, if present, by three nucleotides
    at the threeprime end. Use in cases when annotation files exclude the stop
    codon from the annotated CDS.
    
    Parameters
    ----------
    tx : |Transcript|
        query transcript
        
    Returns
    -------
    |Transcript|
        |Transcript| with same attributes as `tx`, but with CDS extended by one codon
    
    Raises
    ------
    IndexError
        if a three prime UTR is defined that terminates before the complete stop codon
    """
    cdef:
        long start, end, new_start, new_end
        int n
        GenomicSegment span, new_seg, old_seg
        list segments
        dict attr
        Transcript new_tx

    if tx.cds_genome_end is not None and tx.cds_genome_start is not None:
        new_start = start = <long>tx.cds_genome_start
        new_end = end = <long>tx.cds_genome_end
        segments  = copy.copy(tx._segments)
        attr = copy.deepcopy(tx.attr)
        span = tx.spanning_segment
        chrom = span.chrom
        if span.c_strand == forward_strand:
            # if cds at end of transcript, both need to grow by three
            if span.end == end:
                new_end = end + 3
                old_seg = segments[-1]
                new_seg = GenomicSegment(chrom,old_seg.start,new_end,span.strand)
                segments[-1] = new_seg
            else:
                # if new end will be end of transcript, set manually
                # because key won't be in position hash
                if tx.cds_end + 3 == tx.length:
                    new_end = span.end
                else:
                    # position is internal, so we fetch
                    new_end = tx.c_get_genomic_coordinate(tx.cds_end+3,True)
                    # correct for rare cases in which stop codons and at the end
                    # of a given exon interval, and are falsely mapped to first position
                    # of next exon. Instead, we should map to 1 + the genomic position
                    # of the end of the given exon.
                    for n,segment in enumerate(tx):
                        if new_end == segment.start:
                            new_end = tx[n-1].end
        else:
            # if cds starts at beginning of transcript, both need to grow by three
            if span.start == start:
                new_start = start - 3
                old_seg = segments[0]
                new_seg = GenomicSegment(chrom,new_start,old_seg.end,span.strand)
                segments[0] = new_seg
            else:
                # otherwise, fetch
                # minus one here might look funny, but is due to minus-strand
                # feature-ness
                new_start = tx.c_get_genomic_coordinate(tx.cds_end + 3 - 1,True)

        # since `segments` come from an existing transcript, we know they
        # non-overlapping. So, we can bypass checks in __cinit__ and go through
        # faster set_segments()
        new_tx = Transcript()
        new_tx._set_segments(segments)
        new_tx.attr = attr
        new_tx.cds_genome_start = new_start
        new_tx.cds_genome_end = new_end
        new_tx._update_cds()
        return new_tx
    else:
        return tx


#==============================================================================
# Import/export from BED (used by readers)
#==============================================================================
    


cdef dict get_attr_from_bed_column_tuples(list items, int num_bed_columns, int total_columns, list column_tuples):
    """Get SegmentChain attributes from a `BED`_ line with extra columns specified as a list of (column_name,formatter_func)"""
    cdef:
        dict         attr   = get_standard_bed_attr(items, num_bed_columns)
        int          i
        list         xorder = []
    
    for i in range(num_bed_columns,total_columns):
        name, formatter = column_tuples[i-num_bed_columns]
        xorder.append(name)
        try: 
            attr[name] = formatter(items[i])
        except:
            warn("SegmentChain.from_bed(): Could not format '%s' as '%s'. Defaulting to type 'str'." % (items[i],str(formatter.__name__)),DataWarning)
            attr[name] = items[i]

    attr["_bedx_column_order"] = xorder
    return attr

cdef dict get_attr_from_bed_column_names(list items, int num_bed_columns, list column_names): 
    """Get SegmentChain attributes from a `BED`_ line with extra column names given by a list of column names"""
    cdef:
        dict         attr   = get_standard_bed_attr(items, num_bed_columns)
        int          n
        str          name
        
    for n, name in enumerate(column_names):
        attr[name] = items[n+num_bed_columns]

    attr["_bedx_column_order"] = column_names
    return attr

cdef dict get_attr_from_bed_column_number(list items, int num_bed_columns, int num_extra_columns): 
    """Get SegmentChain attributes from a `BED`_ line, given a number of extra columns"""
    cdef:
        dict         attr   = get_standard_bed_attr(items, num_bed_columns)
        int          i
        str          name
        list         xorder = []
    
    for i in range(num_extra_columns):
        name = "custom%s" % i
        attr[name] = items[i+num_bed_columns]
        xorder.append(name)

    attr["_bedx_column_order"] = xorder
    return attr

cdef dict get_standard_bed_attr(list items, int num_bed_columns):
    """Get SegmentChain attributes from standard BED columns in `BED`_ line"""
    cdef:
        str     chrom, strand, default_id
        long    chrom_start, thickstart, thickend
        str     KEY
        object  DEFAULT, _
        tuple   tup
        int     i
        dict    attr, _default_bed_columns

    if num_bed_columns < 3:
        raise ValueError("BED format requires at least 3 columns. Found only %s.\n\t    %s" % (num_bed_columns,items))
    
    chrom         = items[0]
    strand        = "." if num_bed_columns < 6 else items[5]
    chrom_start   = long(items[1])
    chrom_end     = long(items[2])
    
    default_id  = "%s:%s-%s(%s)" % (chrom,chrom_start,chrom_end,strand)

    _default_bed_columns = {
        3 :  ("ID",         default_id,                  str),
        4 :  ("score",      numpy.nan,                   float),
        #5 :  ("strand",    ".", strand),
        6 :  ("thickstart", -1,                          long),
        7 :  ("thickend",   -1,                          long),
        8 :  ("color",      "0,0,0",                     str),
        9 :  ("blocks",     "1",                         int),
        10 : ("blocksizes", str(chrom_end - chrom_start),str),
        11 : ("blockstarts","0",                         str),
    }

    # dict mapping optional bed column to tuple of (Name,default value)
    # these values are used if any optional columns 4-12 are ommited
    # set attr defaults in case we're dealing with BED4-BED9 format
    attr = { KEY : DEFAULT for KEY,DEFAULT,_ in _default_bed_columns.values() }

    for i in (3,4,6,7,8,9,10,11):
        # populate attr with real values from BED columns that are present
        if i >= num_bed_columns:
            break
        
        tup = _default_bed_columns[i]
        key     = tup[0]
        default = tup[1]
        func    = tup[2]
        try:
            attr[key] = func(items[i])
        except IndexError:
            # fall back to default value if real value isn't present
            attr[key] = default
        except ValueError:
            # fall back to default value if real value cannot be parsed
            warn("get_standard_bed_attr: Could not format column %s with '%s'. Falling back to default value '%s'." % (items[i],str(func.__name__),default),DataWarning)
            attr[key] = default

    # convert color to hex string
    try:
        attr["color"] = get_str_from_rgb255(tuple([int(X) for X in attr["color"].split(",")]))
    except ValueError:
        attr["color"] = "#000000"

    # sanity check on thickstart and thickend
    thickstart = attr["thickstart"]
    thickend   = attr["thickend"]

    if thickstart == thickend or None in (thickstart, thickend) \
    or thickstart < 0 or thickend < 0:
        attr["thickstart"] = attr["thickend"] = chrom_start
    
    # stash these temporarily. I don't like this, but hey
    attr["chrom"]       = chrom
    attr["chrom_start"] = chrom_start
    attr["strand"]      = strand
    return attr

cdef dict get_attr_from_bed(str line,object extra_columns=0):
    """See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    for more details.

    Parameters
    ----------
    line
        Line from a `BED`_ file, containing 4 or more columns

    extra_columns: int or list optional
        Extra, non-BED columns in :term:`Extended BED`_ format file corresponding
        to feature attributes. This is common in `ENCODE`_-specific `BED`_ variants.
        
        if `extra-columns` is:
        
          - an :class:`int`: it is taken to be the
            number of attribute columns. Attributes will be stored in
            the `attr` dictionary of the |SegmentChain|, under names like
            `custom0`, `custom1`, ... , `customN`.

          - a :class:`list` of :class:`str`, it is taken to be the names
            of the attribute columns, in order, from left to right in the file.
            In this case, attributes in extra columns will be stored under
            their respective names in the `attr` dict.

          - a :class:`list` of :class:`tuple`, each tuple is taken
            to be a pair of `(attribute_name, formatter_func)`. In this case,
            the value of `attribute_name` in the `attr` dict of the |SegmentChain|
            will be set to `formatter_func(column_value)`.
        
        (Default: 0)

    Returns
    -------
    dict
        Dictionary of attributes
    """
    cdef:
        list  items         = line.strip("\n").split("\t")
        int   total_columns = len(items)
        int   num_bed_columns, num_extra_columns
        dict  attr
        set   types
    
    if extra_columns == 0:
        attr = get_standard_bed_attr(items, len(items))
        
    elif isinstance(extra_columns,int):
        attr = get_attr_from_bed_column_number(items, total_columns - extra_columns, extra_columns)
        
    elif isinstance(extra_columns,(list,tuple)):
        if len(extra_columns) == 0:
            attr = get_standard_bed_attr(items, len(items))
        else:
            types = set([type(X) for X in extra_columns])
            num_extra_columns = len(extra_columns)
            num_bed_columns   = total_columns - num_extra_columns
            
            if len(types) > 1:
                raise ValueError("get_attr_from_bed: List of `extra_columns` contains mixed types. All items need to be the same type, either str or tuple.")
            elif types == { str }:
                attr = get_attr_from_bed_column_names(items, num_bed_columns, extra_columns)
            elif types == { tuple }:
                if all([len(X) == 2 for X in extra_columns]) == False:
                    raise ValueError("get_attr_from_bed: Cannot make SegmentChain from BED input: if a list, extra_columns must be a list of tuples of (column_name,formatter_func)")
                attr = get_attr_from_bed_column_tuples(items, num_bed_columns, total_columns, extra_columns)
            else:
                raise TypeError("get_attr_from_bed: Couldn't figure out types from input '%s' of types '%s'" % (extra_columns,types))
    else:
        raise TypeError("get_attr_from_bed: Cannot make SegmentChain from BED input: extra_columns must be an int or list. Got a %s" % type(extra_columns))

    if attr is None:
        raise ValueError("Could not parse BED line:\n\t    '%s'" % line)
    
    return attr

cdef list get_segments_from_bed(dict attr):
    """Construct a list of |GenomicSegments| given a dictionary of attributes,
    as made from :meth:`get_attr_from_bed` or other methods above. i.e.
    `attr` must contain the following keys: `chrom`, `strand`, `chrom_start`,
    `blocksizes`, `blockstarts`, `blocks`
    """
    cdef:
        str  chrom       = attr["chrom"]
        str  strand      = attr["strand"]
        list segments    = []
        long chrom_start = attr["chrom_start"]
        
    # convert blocks to GenomicSegments
    num_frags   = int(attr["blocks"])
    seg_sizes   = attr["blocksizes"].strip(",").split(",")
    seg_offsets = attr["blockstarts"].strip(",").split(",")
    for i in range(0,num_frags):
        seg_start = chrom_start + int(seg_offsets[i])
        seg_end   = seg_start   + int(seg_sizes[i])
        segments.append(GenomicSegment(chrom,seg_start,seg_end,strand))

    # clean up attr
    for key in ("blocks","blocksizes","blockstarts","chrom_start","strand","chrom"):
        attr.pop(key)
    
    return segments


#==============================================================================
# Helpers
#==============================================================================

cdef void nonecheck(object obj,str place, str valname):
    """Propagate errors if incoming objects are None"""
    if obj is None:
        raise ValueError("%s: value of %s cannot be None" % (place,valname))


cdef bint check_segments(SegmentChain chain, tuple segments) except False:
    cdef:
        GenomicSegment seg, seg0
        GenomicSegment span = chain.spanning_segment
        str my_chrom  = span.chrom
        Strand my_strand = span.c_strand
        int i = 0
        int length = len(segments)
        bint prob = False

    if length > 0:
        msg = "SegmentChain.add_segments: incoming segments (%s) mismatch chain '%s'"
        seg0 = segments[0]
        if len(chain._segments) == 0:
            my_chrom  = seg0.chrom
            my_strand = seg0.c_strand

        while i < length:
            seg = segments[i]
            if seg.chrom != my_chrom:
                msg += "; wrong and/or multiple chromosomes"
                prob = True
                break
                
            if seg.c_strand != my_strand:
                msg += "; wrong and/or multiple strands"
                prob = True
                break

            i += 1

        if prob == True: 
            msg = msg % (", ".join([str(X) for X in segments]),chain)
            raise ValueError(msg)

    return True

cdef ExBool chain_richcmp(SegmentChain chain1, SegmentChain chain2, int cmpval) except bool_exception:
    """Helper method for rich comparisons (==, !=, <, >, <=, >=) of |SegmentChains|

    Parameters
    ----------
    chain1 : SegmentChain
        First chain

    chain2 : SegmentChain
        Second chain

    cmpval : int
        Integer code for comparison (from Cython specification for ``__richcmp__()`` methods:

        ========  ==============
        **Code**  **Comparison**
        --------  --------------
        0         <
        1         <=
        2         ==
        3         !=
        4         >
        5         >=
        ========  ==============
 

    Returns
    -------
    ExBool
        true if `chain1 'comparison' chain2` is `True`, otherwise false
    """
    cdef:
        sspan = chain1.spanning_segment
        ospan = chain2.spanning_segment
        str sname, oname
        long slen, olen

    if cmpval == EQ:
        if len(chain1) == 0 or len(chain2) == 0:
            return false
        else:
            if sspan.chrom == ospan.chrom and\
               sspan.c_strand == ospan.c_strand and\
               chain1._segments == chain2._segments:
                    return true
            return false
    elif cmpval == NEQ:
        return true if chain_richcmp(chain1,chain2,EQ) == false else true
    elif cmpval == LT:
        if sspan < ospan:
            return true
        elif sspan > ospan:
            return false
        else:
            slen = chain1.length
            olen = chain2.length
            if slen < olen:
                return true
            elif slen > olen:
                return false
            else:
                sname = chain1.get_name()
                oname = chain2.get_name()
                if sname < oname:
                    return true
                else:
                    return false
    elif cmpval == LEQ:
        return true if chain_richcmp(chain1,chain2,LT) == true \
                    or chain_richcmp(chain1,chain2,EQ) == true \
                    else false
    elif cmpval == GT:
        return chain_richcmp(chain2,chain1,LT)
    elif cmpval == GEQ:
        return true if chain_richcmp(chain2,chain1,LT) == true \
                    or chain_richcmp(chain1,chain2,EQ) == true \
                    else false
    else:
        raise ValueError("This operation is not defined for SegmentChains.")

cdef ExBool transcript_richcmp(Transcript chain1, Transcript chain2, int cmpval) except bool_exception:
    """Helper method for rich comparisons (==, !=, <, >, <=, >=) of |Transcripts|

    Parameters
    ----------
    chain1 : Transcript
        First chain

    chain2 : Transcript
        Second chain

    cmpval : int
        Integer code for comparison (from Cython specification for ``__richcmp__()`` methods:

        ========  ==============
        **Code**  **Comparison**
        --------  --------------
        0         <
        1         <=
        2         ==
        3         !=
        4         >
        5         >=
        ========  ==============
 

    Returns
    -------
    ExBool
        true if `chain1 'comparison' chain2` is `True`, otherwise false
    """
    cdef:
        sspan = chain1.spanning_segment
        ospan = chain2.spanning_segment

    if cmpval == EQ:
        if chain_richcmp(chain1,chain2,cmpval) == true and \
               (
                   (chain1.cds_genome_start is None and chain2.cds_genome_start is None and\
                    chain1.cds_genome_end is None and chain2.cds_genome_end is None) or\
                   (chain1.cds_genome_start == chain2.cds_genome_start and\
                    chain1.cds_genome_end == chain2.cds_genome_end)
                ):
                    return true
        return false
    elif cmpval == NEQ:
        return true if transcript_richcmp(chain1,chain2,EQ) == false else true
    elif cmpval == LT:
        return chain_richcmp(chain1,chain2,cmpval)
    elif cmpval == LEQ:
        return true if chain_richcmp(chain1,chain2,LT) == true \
                    or transcript_richcmp(chain1,chain2,EQ) == true \
                    else false
    elif cmpval == GT:
        return chain_richcmp(chain2,chain1,LT)
    elif cmpval == GEQ:
        return true if chain_richcmp(chain2,chain1,LT) == true \
                    or transcript_richcmp(chain1,chain2,EQ) == true\
                    else false
    else:
        raise ValueError("This operation is not defined for Transcripts.")



#==============================================================================
# Classes
#==============================================================================

cdef class GenomicSegment:
    """
    GenomicSegment(chrom, start, end, strand)
    
    Building block for |SegmentChain|: a continuous segment of the genome
    defined by a chromosome name, start coordinate,  end coordinate, and strand.

    Attributes
    ----------
    chrom : str
        Name of chromosome

    start : int
        0-indexed, left most position of segment

    end : int
        0-indexed, half-open right most position of segment

    strand : str
        Chromosome strand (`'+'`, `'-'`, or `'.'`)
   

    Examples
    --------

    |GenomicSegments| sort lexically by chromosome, start position, end position,
    and finally strand::
        
        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrB",0,10,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",75,100,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",55,75,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",50,150,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") < GenomicSegment("chrA",50,100,"-")
        True


    They also provide a few convenience methods for containment or overlap. To be
    contained, a segment must be on the same chromosome and strand as its container,
    and its coordinates must be within or equal to its endpoints::

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",25,100,"+")
        True

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",50,100,"+")
        True
       
        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",25,100,"-")
        False

        >>> GenomicSegment("chrA",50,100,"+") in GenomicSegment("chrA",75,200,"+")
        False

    Similarly, to overlap, |GenomicSegments| must be on the same strand and 
    chromosome.
    """
    def __cinit__(self, str chrom, long start, long end, str strand):
        """Create a |GenomicSegment|
        
        Parameters
        ----------
        chrom : str
            Chromosome name
        
        start : int
            0-indexed, leftmost coordinate of feature
        
        end : int
            0-indexed, half-open rightmost coordinate of feature
            Must be >= `start`
        
        strand : str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
        """
        if end < start:
            raise ValueError("GenomicSegment: start coordinate (%s) must be >= end (%s)." % (start,end))
        my_chrom = chrom
        self.chrom  = my_chrom
        self.start  = start
        self.end    = end
        self.c_strand = str_to_strand(strand)
    
    def __reduce__(self): # enables pickling
        return (GenomicSegment,(self.chrom,self.start,self.end,self.strand))

    def __repr__(self):
        return "<%s %s:%s-%s strand='%s'>" % ("GenomicSegment",
                                              self.chrom,
                                              self.start,
                                              self.end,
                                              self.strand)
    
    def __str__(self):
        return "%s:%s-%s(%s)" % (self.chrom, self.start, self.end, self.strand)
   
    def __hash__(self):
        # generate unique hash from a tuple, which are hashable
        # and already implement hashing for us
        return hash((self.__class__.__name__,self.chrom,self.start,self.end,self.strand))

    @staticmethod
    def from_str(str inp):
        """Construct a |GenomicSegment| from its ``str()`` representation
        
        Parameters
        ----------
        inp : str
            String representation of |GenomicSegment| as `chrom:start-end(strand)`
            where `start` and `end` are in 0-indexed, half-open coordinates
        
        Returns
        -------
        |GenomicSegment|
        """
        cdef:
            long start, end
            str chrom, strand, s_start, s_end
        chrom, s_start, s_end, strand = segpat.search(inp).groups()
        start = long(s_start)
        end   = long(s_end) 
        return GenomicSegment(chrom,start,end,strand)
    
    def __copy__(self):
        cdef str newchrom, newstrand
        cdef long newstart, newend
        newchrom = copy.copy(self.chrom)
        newstrand = self.strand
        newstart = self.start
        newend = self.end
        return GenomicSegment(newchrom,newstart,newend,newstrand)

    def __deepcopy__(self,memo):
        return self.__copy__()

    def __len__(self):
        """Return length, in nucleotides, of |GenomicSegment|"""
        cdef long len
        len = self.end - self.start
        return len

    def __richcmp__(self,GenomicSegment other, int cmptype):
        if other is None or not isinstance(other,GenomicSegment):
            return False

        return self._cmp_helper(other,cmptype)

    cpdef bint _cmp_helper(self,GenomicSegment other,int cmptype):
        nonecheck(other,"GenomicSegment eq/neq","other")
        schrom = self.chrom
        ochrom = other.chrom
        if cmptype == EQ:
            return self.start     == other.start and\
                   self.end       == other.end and\
                   schrom         == ochrom and\
                   self.c_strand  == other.c_strand
        elif cmptype == NEQ:
            return self._cmp_helper(other,EQ) == False
        elif cmptype == LT:
            if schrom < ochrom:
                return True
            elif schrom > ochrom:
                return False
            elif schrom == ochrom:
                sstart = self.start
                ostart = other.start
                if sstart < ostart:
                    return True
                elif sstart > ostart:
                    return False
                elif sstart == ostart:
                    send = self.end
                    oend = other.end
                    if send < oend:
                        return True
                    elif send > oend:
                        return False
                    elif send == oend:
                        return self.c_strand < other.c_strand
        elif cmptype == GT:
            return other._cmp_helper(self,LT) # (other < self)
        elif cmptype == LEQ:
            return self._cmp_helper(other,LT) or self._cmp_helper(other,EQ)
        elif cmptype == GEQ:
            return other._cmp_helper(self,LT) or self._cmp_helper(other,EQ)
        else:
            raise AttributeError("Comparison operation not defined")
    
    def __contains__(self,other):
        """Test whether this segment contains `other`, where containment is
        defined as all positions in `other` being present in self, when both
        `self` and `other` share the same chromosome and strand.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        nonecheck(other,"GenomicSegment.__eq__","other")
        return self.contains(other)

    cpdef bint contains(self,GenomicSegment other):
        """Test whether this segment contains `other`, where containment is
        defined as all positions in `other` being present in self, when both
        `self` and `other` share the same chromosome and strand.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        nonecheck(other,"GenomicSegment.contains","other")
        return self.chrom == other.chrom and\
               self.c_strand == other.c_strand and\
               (other.start >= self.start and other.end <= self.end and other.end >= other.start)
         
    cpdef bint overlaps(self,GenomicSegment other):
        """Test whether this segment overlaps `other`, where overlap is defined
        as sharing: a chromosome, a strand, and a subset of coordinates.
           
        Parameters
        ----------
        other : |GenomicSegment|
            Query segment
        
        Returns
        -------
        bool
        """
        nonecheck(other,"GenomicSegment.overlaps","other")
        if self.chrom == other.chrom and self.c_strand == other.c_strand:
            if (self.start >= other.start and self.start < other.end) or\
               (other.start >= self.start and other.start < self.end):
                    return True

        return False

    cpdef str as_igv_str(self):
        """Format as an IGV location string"""
        return "%s:%s-%s" % (self.chrom, self.start+1, self.end+1)
    
    @staticmethod
    def from_igv_str(str loc_str, str strand="."):
        """Construct |GenomicSegment| from IGV location string
        
        Parameters
        ----------
        igvloc : str
            IGV location string, in format `'chromosome:start-end'`,
            where `start` and `end` are 1-indexed and half-open
            
        strand : str
            The chromosome strand (`'+'`, `'-'`, or `'.'`)
            
        Returns
        -------
        |GenomicSegment|
        """
        cdef:
            str chrom, s_start, s_end
            long start, end

        chrom,s_start,s_end = igvpat.search(loc_str).groups()
        start = long(s_start) - 1
        end   = long(s_end) - 1
        return GenomicSegment(chrom,start,end,strand)
    
    property start:
        """Zero-indexed (Pythonic) start coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.start

    property end:
        """Zero-indexed, half-open (Pythonic) end coordinate of |GenomicSegment|"""
        def __get__(self):
            return self.end

    property chrom:
        """Chromosome where |GenomicSegment| resides"""
        def __get__(self):
            return self.chrom

    property strand:
        """Strand of |GenomicSegment|

          - '+' for forward / Watson strand
          - '-' for reverse / Crick strand
          - '.' for unstranded / both strands
        """
        def __get__(self):
            return strand_to_str(self.c_strand)

    property c_strand:
        def __get__(self):
            return self.c_strand

cdef class SegmentChain(object):
    """
    SegmentChain(*segments,**attributes)
    
    Base class for genomic features, composed of zero or more |GenomicSegments|.
    |SegmentChains| can therefore model discontinuous, features -- such as
    multi-exon transcripts or gapped alignments -- in addition, to continuous
    features.
    
    Numerous convenience functions are supplied for:
    
      - converting between coordinates relative to the genome and relative
        to the internal coordinates of a spliced |SegmentChain|
        
      - fetching genomic sequence, read alignments, or count data, accounting
        for splicing of the segments, and, in the case of reverse-strand features,
        reverse-complementing
      
      - slicing or fetching sub-regions of a |SegmentChain|
      
      - testing equality, inequality, overlap, containment, coverage of, or
        sharing of segments with other |SegmentChain| or |GenomicSegment| objects
    
      - import/export to `BED`_, `PSL`_, `GTF2`_, and `GFF3`_ formats,
        for use in other software packages or in a genome browser.
    
    Intervals are sorted from lowest to greatest starting coordinate on their
    reference sequence, regardless of strand. Iteration over the SegmentChain
    will yield intervals from left-to-right in the genome.


    Parameters
    ----------
    *segments : |GenomicSegment|
        0 or more GenomicSegments on the same strand
    
    **attr : keyword arguments
        Arbitrary attributes, including, for example:
    
        ====================    ============================================
        **Attribute**           **Description**
        ====================    ============================================
        ``type``                A feature type used for `GTF2`_/`GFF3`_ export
                                of each interval in the |SegmentChain|. (Default: `'exon'`)
        
        ``ID``                  A unique ID for the |SegmentChain|.

        ``transcript_id``       A transcript ID used for `GTF2`_ export

        ``gene_id``             A gene ID used for `GTF2`_ export
        ====================    ============================================    


    Attributes
    ----------
    spanning_segment : |GenomicSegment|
        A |GenomicSegment| spanning the endpoints of the |SegmentChain|

    strand : str
        The chromosome strand (`'+'`, `'-'`, or `'.'`)

    chrom : str
        Name of the chromosome on which the |SegmentChain| resides

    attr : dict
        Any miscellaneous attributes or annotation data
        
    segments : list
        List of |GenomicSegments| in the chain
    
    mask_segments : list
        List of |GenomicSegments| that mask positions in the chain.

    See Also
    --------
	Transcript
        Transcript subclass, additionally providing richer `GTF2`_, `GFF3`_,
        and `BED`_ export, as well as methods for fetching coding regions
        and UTRs as subsegments
    """

    def __cinit__(self,*segments,**attr):
        """Create an |SegmentChain| from zero or more |GenomicSegment| objects
        
        Example::
        
            >>> seg1 = GenomicSegment("chrI",2000,2500,"+")
            >>> seg2 = GenomicSegment("chrI",10000,11000,"+")
            >>> chain = SegmentChain(seg1,seg2,ID="example_chain",type="mRNA")
            
        
        Parameters
        ----------
        *segments : |GenomicSegment|
            0 or more GenomicSegments on the same strand
        
        **attr : dict
            Arbitrary attributes, including, for example:
        
            ====================    ============================================
            **Attribute**           **Description**
            ====================    ============================================
            ``type``                A feature type used for `GTF2`_/`GFF3`_ export
                                    of each interval in the |SegmentChain|. (Default: `'exon'`)
            
            ``ID``                  A unique ID for the |SegmentChain|.

            ``transcript_id``       A transcript ID used for `GTF2`_ export

            ``gene_id``             A gene ID used for `GTF2`_ export
            ====================    ============================================
        """
        cdef:
            int num_segs = len(segments)
            list new_segments

        if "type" not in attr:
            attr["type"] = "exon"

        self.attr = attr
        self._mask_segments = None 
        self._segments      = []
        self.length         = 0
        self.masked_length  = 0
        self.spanning_segment = NullSegment
        self._position_mask = None
        self._position_hash = None
        if num_segs == 0:
            pass
        elif num_segs == 1:
            self._set_segments(list(segments))
        else:
            check_segments(self,segments)
            new_segments = merge_segments(list(segments))
            self._set_segments(new_segments)

    cdef bint _set_segments(self, list segments) except False:
        """Set `self._segments` assuming `segments` have passed the following criteria:

            1. They don't overlap

            2. They all have the same chromosome and strand

            3. If `len(self)` > 0, all segments have the same chromosome
               and strand as `self`
            
            4. They are sorted

        Ordinarily, these criteria are checked by :func:`~plastid.genomics.c_segmentchain.check_segments`
        and either :meth:`~plastid.genomics.c_segmentchain.SegmentChain.__cinit__`
        or :meth:`~plastid.genomics.c_segmentchain.SegmentChain.add_segments`.

        Occasionally (e.g. in the case of unpickling) it is safe to bypass these
        checks.

        The following properties are updated:

            `self.length`
                Length in nucleotides of `self`

            `self._segments`
                |GenomicSegments| in `self`

            `self.spanning_segment`
                |GenomicSegment| spanning entire length of `self`

        And :meth:`SegmentChain.c_reset_masks` is called to invalidate/remove
        any masks.


        Parameters
        ----------
        list
            List of |GenomicSegments|
        """
        cdef:
            int num_segs = len(segments)
            GenomicSegment segment, seg0
            long length = sum([len(X) for X in segments])
        
        self._segments = segments
        self.length = length
        self._position_hash = None
        self.c_reset_masks()
        
        if num_segs == 0:
            self.spanning_segment = NullSegment
        elif num_segs == 1:
            self.spanning_segment = self._segments[0]
        else:
            segs = self._segments
            seg0 = segs[0]
            self.spanning_segment = GenomicSegment(seg0.chrom,
                                                   seg0.start,
                                                   segs[-1].end,
                                                   seg0.strand)

        return True

    cdef array.array _get_position_hash(self):
        """Populate hash mapping coordinates in `self` to corresponding genomic coordinates,
        leftmost first (i.e. position 0 is that with the lowest corresponding genomic coordinate, 
        regardless of `self.strand`).

        Results are saved in `self._position_hash`. If `self._position_hash` is not `None`,
        it will be returned instead of being repopulated. To force repopulation, first set
        `self._position_hash` to `None`.

        Returns
        -------
        :class:`array.array`
            Array of `length(self)`, mapping coordinates in `self` to genomic coordinates,
            leftmost first (i.e. in genomic order, not transcript order).
        """
        if self._position_hash is not None:
            return self._position_hash

        cdef:
            GenomicSegment segment
            long x, c
            array.array my_hash = array.clone(hash_template,self.length,False)  
            long [:] my_view = my_hash
        
        self._position_hash = my_hash

        c = 0
        for segment in self._segments:
            x = segment.start
            while x < segment.end:
                my_view[c] = x
                c += 1
                x += 1            

        return self._position_hash

    def sort(self): # this should never need to be called, now that _segments and _mask_segments are managed
        self._segments.sort()
        if self._mask_segments is not None:
            self._mask_segments.sort()
        
    property chrom:
        """Chromosome the SegmentChain resides on"""
        def __get__(self):
            return self.spanning_segment.chrom

    property strand:
        """Strand of the SegmentChain"""
        def __get__(self):
            return strand_to_str(self.spanning_segment.c_strand)

    property c_strand:
        def __get__(self):
            return self.spanning_segment.c_strand

    property segments:
        """Copy of list of |GenomicSegments| that comprise `self`.
        Changing this list will do nothing to `self`.
        """
        def __get__(self):
            return copy.deepcopy(self._segments)

    property mask_segments:
        """Copy of list of |GenomicSegments| representing regions masked in `self.`
        Changing this list will do nothing to the masks in `self`.
        """
        def __get__(self):
            if self._mask_segments is None:
                return []

            return copy.deepcopy(self._mask_segments)

    def __reduce__(self): # enable pickling
        return (self.__class__,tuple(),self.__getstate__())

    def __getstate__(self): # save state for pickling
        cdef:
            list segstrs  = [str(X) for X in self._segments]
            list maskstrs = []
            dict attr = self.attr

        if self._mask_segments is not None:
            maskstrs = [str(X) for X in self._mask_segments]


        return (segstrs,maskstrs,attr)

    def __setstate__(self,state): # revive state from pickling
        cdef:
            list segs  = [GenomicSegment.from_str(X) for X in state[0]]
            list masks = [GenomicSegment.from_str(X) for X in state[1]]
            dict attr = state[2]

        self.attr = attr
        self._set_segments(segs)
        self._set_masks(masks)

    def __copy__(self):
        chain2 = SegmentChain()
        chain2.attr = self.attr
        chain2._set_segments(self._segments)
        if self._mask_segments is not None:
            chain2._set_masks(self._mask_segments)
        return chain2

    def __deepcopy__(self,memo):
        chain2 = SegmentChain()
        chain2.attr = copy.deepcopy(self.attr)
        chain2._set_segments(self._segments)
        if self._mask_segments is not None:
            chain2._set_masks(copy.deepcopy(self._mask_segments))
        return chain2

    def __richcmp__(self, object other, int cmpval):
        """Test whether `self` and `other` are equal or unequal. Equality is defined as
        identity of positions, chromosomes, and strands. Two |SegmentChain| with
        zero intervals, by convention, are not equal.
           
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
        """
        if isinstance(other,SegmentChain):
            return chain_richcmp(self,other,cmpval) == true
        elif isinstance(other,GenomicSegment):
            return chain_richcmp(self,SegmentChain(other),cmpval) == true 
        else:
            raise TypeError("SegmentChain eq/ineq/et c is only defined for other SegmentChains or GenomicSegments.")

    def __repr__(self):
        cdef:
            str sout = "<%s segments=%s" % (self.__class__.__name__, len(self._segments))
            GenomicSegment span = self.spanning_segment

        if len(self) > 0:
            sout += " bounds=%s:%s-%s(%s)" % (span.chrom,
                                              span.start,
                                              span.end,
                                              span.strand)
            sout += " name=%s" % self.get_name()
        sout += ">"
        return sout

    def __str__(self):
        """String representation of |SegmentChain|. Inverse of :py:meth:`.from_str`
        Chains are represented as:
        
            `'chrom_name:segment1_start-segment1_end^segment2_start-segment2_end^...^(strand)'`
        
        Where all coordinates are zero-indexed and half-open. Masked segments
        are not saved in this representation; nor are attributes in `self.attr`
        """ 
        if len(self) > 0:
            ltmp = ["%s-%s" % (segment.start, segment.end) for segment in self]
            stmp = "^".join(ltmp)
            sout = "%s:%s(%s)" % (self.chrom,stmp,self.strand)
        else:
            sout = "na"
        return sout

    def __getitem__(self,int index):
        """Fetch a |GenomicSegment| from the |SegmentChain|
        
        Parameters
        ----------
        index : int
            Index of interval to select, from left-to-right in genome
        
        Returns
        -------
        |GenomicSegment|
        """
        return self._segments[index]
           
    def __iter__(self):
        """Interation over each |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return iter(self._segments)
    
    def __next__(self):
        """Return next |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return next(self._segments)
    
    def next(self):
        """Return next |GenomicSegment| in the |SegmentChain|,
        from left to right on the chromsome"""
        return self.__next__()
    
    def __len__(self):
        """Return the number of |GenomicSegments| in the |SegmentChain|"""
        return len(self._segments)
    
    def __contains__(self, object other):
        """Tests whether |SegmentChain| contains another |SegmentChain|
        or |GenomicSegment|. Containment is defined for each type as follows:
        
        =====================   =======================================================
        **Type of `other`**     **True if**
        ---------------------   -------------------------------------------------------
        |SegmentChain|          If `self` and `other` both contain more than one segment,
                                all segment-segment junctions in `other` must be       
                                represented in `self`, in identical order, and all
                                positions covered by `other` must be present in `self`.
                                If `other` contains one segment, it must be fully
                                contained by one segment in `self`.                                  

        |GenomicSegment|        `other` must be completely contained        
                                within one of the segments in `self`                       
        =====================   =======================================================
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        
        Returns
        -------
        bool
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.__contains__(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_contains(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.unstranded_overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("The 'in'/containment operator is only defined for GenomicSegments and SegmentChains")
     
    cdef ExBool c_contains(self, SegmentChain other) except bool_exception:
        cdef:
            set selfpos, opos
            list myjuncs, ojuncs
            int i, mystart
            bint found
            GenomicSegment sspan = self.spanning_segment
            GenomicSegment ospan = other.spanning_segment
            list mine = self._segments

        if len(self) == 0 or len(other) == 0:
            return false
        elif sspan.chrom != ospan.chrom:
            return false
        elif sspan.c_strand != ospan.c_strand:
            return false
        elif other.length > self.length:
            return false

        # TODO: better algorithm here
        elif len(other) == 1:
            for segment in self:
                if segment.contains(other[0]):
                    return true
            return false 
        else:
            selfpos = self.c_get_position_set()
            opos    = other.c_get_position_set()
            
            if opos & selfpos == opos:
                myjuncs = self.get_junctions()
                ojuncs  = other.get_junctions()
    
                found = False
                for i, myjunc in enumerate(myjuncs):
                    if ojuncs[0] == myjunc:
                        mystart = i
                        found   = True
                        break
            else:
                return false
            
            if found == True:
                return true if ojuncs == myjuncs[mystart:mystart+len(ojuncs)] else false
            else:
                return false

    def shares_segments_with(self, object other):
        """Returns a list of |GenomicSegment| that are shared between `self` and `other`
           
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        list
            List of |GenomicSegments| common to `self` and `other`
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        if isinstance(other,GenomicSegment):
            return self.c_shares_segments_with(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            return self.c_shares_segments_with(other)
        else:
            raise TypeError("SegmentChain.shares_segments_with() is defined only for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef list c_shares_segments_with(self, SegmentChain other):
        cdef:
            list shared
            GenomicSegment segment
        if self.chrom != other.chrom or self.c_strand != other.c_strand:
            return []
        else:
            shared = []
            for segment in other:
                if segment in self._segments:
                    shared.append(segment)
            return shared
   
    def  unstranded_overlaps(self, object other):
        """Return `True` if `self` and `other` share genomic positions
        on the same chromosome, regardless of their strands
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome, False otherwise. Strands of `self` and `other` need
            not match
            
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,SegmentChain):
            return self.c_unstranded_overlaps(other) == true
        elif isinstance(other,GenomicSegment):
            return self.c_unstranded_overlaps(SegmentChain(other)) == true
        else:
            raise TypeError("SegmentChain.unstranded_overlaps() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))

    cdef ExBool c_unstranded_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            GenomicSegment left
            GenomicSegment right
            str chrom = self.chrom
            list segments
            int i = 0
            int length

        if self.chrom != other.chrom:
            return false

        segments = sorted(self._segments + other._segments)
        length = len(segments)
        if len(segments) > 0:
            left = segments[0]

            while i < length - 1:
                right = segments[i+1]
                if right.start < left.end:
                    return true

                left = right
                i += 1

        return false

    def overlaps(self, object other):
        """Return `True` if `self` and `other` share genomic positions on the same strand
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome and strand; False otherwise.
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.c_overlaps(SegmentChain(other)) == true
        elif isinstance(other,SegmentChain):
            return self.c_overlaps(other) == true
        else:
            raise TypeError("SegmentChain.overlaps() is defined only for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef ExBool c_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            Strand sstrand = self.c_strand
            Strand ostrand = other.c_strand

        if self.c_unstranded_overlaps(other) == true and self.c_strand & other.c_strand > 0:
            return true

        return false

    def antisense_overlaps(self, object other):
        """Returns `True` if `self` and `other` share genomic positions on opposite strands
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share genomic positions on the same
            chromosome but opposite strand; False otherwise.
                    
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,GenomicSegment):
            return self.c_antisense_overlaps(SegmentChain(other))
        elif isinstance(other,SegmentChain):
            retval = self.c_antisense_overlaps(other)
            if retval == bool_exception:
                raise RuntimeError("SegmentChain.antisense_overlaps() errored on chain '%s' and object '%s'" % (self,other))
            return retval == true
        else:
            raise TypeError("SegmentChain.antisense_overlaps() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))

    cdef ExBool c_antisense_overlaps(self, SegmentChain other) except bool_exception:
        cdef:
            Strand sstrand = self.c_strand
            Strand ostrand = other.c_strand
        if self.c_unstranded_overlaps(other) == true:
            if unstranded in (sstrand,ostrand) or sstrand != ostrand:
                return true
        return false

    def covers(self, object other):
        """Return `True` if `self` and `other` share a chromosome and strand,
        and all genomic positions in `other` are present in `self`.
        By convention, zero-length |SegmentChains| are not covered by other
        chains.
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature

         
        Returns
        -------
        bool
            `True` if `self` and `other` share a chromosome and strand, and all
            genomic positions in `other` are present in `self`. Otherwise `False`

        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef ExBool retval
        if isinstance(other,SegmentChain):
            return self.c_covers(other) == true
        elif isinstance(other,GenomicSegment):
            return self.c_covers(SegmentChain(other)) == true
        else:
            raise TypeError("SegmentChain.covers() is only defined for GenomicSegments and SegmentChains. Found %s." % type(other))
 
    cdef ExBool c_covers(self, SegmentChain other) except bool_exception:
        """Return `True` if `self` and `other` share a chromosome and strand,
        and all genomic positions in `other` are present in `self`.
        By convention, zero-length |SegmentChains| are not covered by other
        chains.
        
        
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
         
        Returns
        -------
        bool
            `True` if `self` and `other` share a chromosome and strand, and all
            genomic positions in `other` are present in `self`. Otherwise `False`
        
        Raises
        ------
        TypeError
            if `other` is not a |GenomicSegment| or |SegmentChain|
        """
        cdef:
            GenomicSegment sspan = self.spanning_segment
            GenomicSegment ospan = other.spanning_segment
            list mine = self._segments

        if len(self) == 0 or len(other) == 0:
            return false
        elif sspan.chrom  == ospan.chrom and\
             sspan.c_strand == ospan.c_strand and\
             mine == merge_segments(mine + other._segments):
            return true

        return false
        
#         extra_attr : keyword arguments
#             Values that will be included in the new chain's `attr` dict.
#             These can be used to overwrite values already present.    
    cdef SegmentChain _restrand(self,Strand strand):#,**extra_attr):
        """Create a |SegmentChain| on strand `strand` using coordinates from `self`.
        
        Parameters
        ----------
        strand : Strand
            New strand

                        
        Returns
        -------
        SegmentChain
            |SegmentChain| with new strand
        """
        cdef: 
            str            s_strand     = strand_to_str(strand)
            str            chrom        = self.spanning_segment.chrom
            SegmentChain   new_chain    = SegmentChain()#**extra_attr)
            list           new_segments = []
            GenomicSegment seg
    
        for seg in self._segments:
            new_segments.append(GenomicSegment(chrom,seg.start,seg.end,s_strand))

        new_chain._set_segments(new_segments)
        return new_chain
    
    cpdef SegmentChain get_unstranded(self):#,**extra_attr):
        """Returns an |SegmentChain| antisense to `self`, with empty `attr` dict.

        Returns
        -------
        SegmentChain
            |SegmentChain| antisense to `self`
        """
        return self._restrand(unstranded)#,**extra_attr)
    
    cpdef SegmentChain get_antisense(self):#,**extra_attr):
        """Returns an |SegmentChain| antisense to `self`, with empty `attr` dict.

        Returns
        -------
        SegmentChain
            |SegmentChain| antisense to `self`
        """
        cdef:
            Strand strand = self.spanning_segment.c_strand

        if strand == forward_strand:
            new_strand = reverse_strand
        elif strand == reverse_strand:
            new_strand = forward_strand
        elif strand == unstranded:
            new_strand = unstranded
        else:
            raise ValueError("Strand for SegmentChain '%s' is undefined and has no anti-sense." % self)
        
        return self._restrand(new_strand)#,**extra_attr)

    cdef list c_get_position_list(self):
        """Retrieve a sorted **end-inclusive** numpy array of genomic coordinates in this |SegmentChain|
         
        Returns
        -------
        list
            Genomic coordinates in `self`, as integers, in genomic order
        """
        return self._get_position_hash().tolist()

    def get_position_list(self):
        """Retrieve a sorted **end-inclusive** numpy array of genomic coordinates in this |SegmentChain|
         

        Returns
        -------
        list
            Genomic coordinates in `self`, as integers, in genomic order
        """
        return self._get_position_hash().tolist()

    cdef numpy.ndarray c_get_position_array(self, bint copy):
        """Retrieve a sorted end-inclusive list of genomic coordinates in this |SegmentChain|
        
        Parameters
        ----------
        copy : bool, optional
            If `False` (default), return a view of the |SegmentChain|'s
            internal position mapping. If `True`, return  a copy

        Returns
        -------
        :class:`numpy.ndarray`
            Ggenomic coordinates in `self`, as integers, in genomic order
        """
        cdef numpy.ndarray[LONG_t,ndim=1] positions = numpy.asarray(self._get_position_hash(),dtype=long) 
        if copy == False:
            return positions
        else:
            return copy.deepcopy(positions)

    def get_position_set(self):
        """Retrieve an end-inclusive set of genomic coordinates included in this |SegmentChain|
        
        Returns
        -------
        set
            Set of genomic coordinates, as integers
        """
        return set(self.c_get_position_list())

    cdef set c_get_position_set(self):
        return set(self._get_position_hash().tolist())
   
    cpdef set get_masked_position_set(self):
        """Returns a set of genomic coordinates corresponding to positions in 
        `self` that **HAVE NOT** been masked using :meth:`SegmentChain.add_masks`

        Returns
        -------
        set
            Set of genomic coordinates, as integers
        """
        cdef:
            numpy.ndarray[int,ndim=1] mask 
            numpy.ndarray[LONG_t,ndim=1] positions

        if self._position_mask is None:
            return self.c_get_position_set()
        else:
            mask = numpy.frombuffer(self._position_mask,dtype=numpy.intc)
            positions = numpy.frombuffer(self._get_position_hash(),dtype=LONG)
            return set(positions[mask == 0])
    
    def get_name(self):
        """Returns the name of this |SegmentChain|, first searching through
        `self.attr` for the keys `ID`, `Name`, and `name`. If no value is found
        for any of those keys, a name is generated using :meth:`SegmentChain.__str__`
        
        Returns
        -------
        str
            In order of preference, `ID` from `self.attr`, `Name` from
            `self.attr`, `name` from `self.attr` or ``str(self)`` 
        """
        name = self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self))))
        return name
    
    def get_gene(self):
        """Return name of gene associated with |SegmentChain|, if any, 
        by searching through `self.attr` for the keys `gene_id` and `Parent`.
        If one is not found, a generated gene name for the SegmentChain is 
        made from :py:meth:`get_name`.

        Returns
        -------
        str
            Returns in order of preference, `gene_id` from `self.attr`, 
            `Parent` from `self.attr` or ``'gene_%s' % self.get_name()``
        """
        gene = self.attr.get("gene_id",
               self.attr.get("Parent",
               "gene_%s" % self.get_name()))
        if isinstance(gene,list):
            gene = ",".join(sorted(gene))
            
        return gene
    
    def get_length(self):
        """Return total length, in nucleotides, of `self`
        
        Returns
        -------
        int
        """
        cdef str name = self.__class__.__name__
        warn("%s.get_length() is deprecated and will be removed in future versions. Use %s.length property instead" % (name,name),
             DeprecationWarning)
        return self.length

    def get_masked_length(self):
        """Return the total length, in nucleotides, of positions in `self`
        that have not been masked using :meth:`SegmentChain.add_masks`
        
        Returns
        -------
        int
        """
        cdef str name = self.__class__.__name__
        warn("%s.get_masked_length() is deprecated and will be removed in future versions. Use %s.masked_length property instead" % (name,name),
             DeprecationWarning)
        return self.masked_length

    cdef bint c_add_segments(self, tuple segments) except False:
        """Add |GenomicSegments| to `self`. If there are
        already segments in the chain, the incoming segments must be 
        on the same strand and chromosome as all others present.

        Parameters
        ----------
        segments : tuple
            Tuple of |GenomicSegments|
        """
        cdef:
            int length = len(segments)
            list new_segments

        if length > 0:
            check_segments(self,segments)
            new_segments = list(segments)
            new_segments = merge_segments(new_segments + self._segments)
            self._set_segments(new_segments)

        return True

    def add_segments(self,*segments):
        """Add 1 or more |GenomicSegments| to the |SegmentChain|. If there are
        already segments in the chain, the incoming segments must be 
        on the same strand and chromosome as all others present.

        Parameters
        ----------
        segments : |GenomicSegment|
            One or more |GenomicSegment| to add to |SegmentChain|
        """
        if len(segments) > 0:
            if self._mask_segments is not None and len(self._mask_segments) > 0:
                warn("Segmentchain: adding segments to %s will reset its masks!" % (self),UserWarning)

            self.c_add_segments(segments)
        
    def add_masks(self,*mask_segments):
        """Adds one or more |GenomicSegment| to the collection of masks.
        Masks will be trimmed to the positions of the |SegmentChain|
        during addition.

        Parameters
        ----------
        mask_segments : |GenomicSegment|
            One or more segments, in genomic coordinates, covering positions to
            exclude from return values of :meth:`get_masked_position_set`,
            :meth:`get_masked_counts`, or :meth:`get_masked_length`
		
        See also
        --------
        SegmentChain.get_masks
        SegmentChain.get_masks_as_segmentchain
        SegmentChain.reset_masks
        """
        cdef:
            str my_chrom, my_strand
            GenomicSegment segment
            set positions = set()
            list new_segments, segs

        if len(mask_segments) > 0:
            check_segments(self,mask_segments)
            seg  = mask_segments[0]
            my_chrom  = seg.chrom
            my_strand = seg.strand

            segs = list(mask_segments)
            if self._mask_segments is not None:
                segs += self._mask_segments
               
            # TODO: optimize this- could be way faster
            # add new positions to any existing masks
            for segment in segs:
                positions |= set(range(segment.start,segment.end))
             
            # trim away non-overlapping masks
            positions &= self.c_get_position_set()
            
            # regenerate list of segments from positions, in case some were doubly-listed
            new_segments = positions_to_segments(my_chrom,my_strand,positions)
            self._set_masks(new_segments)
   
    cdef bint _set_masks(self, list segments) except False:
        """Set `self._mask_segments` and update mask hashes, assuming `segments`
        have passed the following criteria:

            1. They don't overlap

            2. They all have the same chromosome and strand as `self`

            3. They are sorted

            4. They are trimmed to the exon boundaries of `self`.

        Ordinarily, these criteria are jointly enforced by
        :func:`~plastid.genomics.c_segmentchain.check_segments`
        and :meth:`~plastid.genomics.c_segmentchain.SegmentChain.add_masks`

        Occasionally (e.g. in the case of unpickling) it is safe (and fast)
        to bypass these checks.

        Parameters
        ----------
        list
            List of |GenomicSegments|
        """
        cdef:
            array.array pmask = array.clone(mask_template,self.length,True)
            int [:] pview = pmask
            list new_segments
            GenomicSegment seg
            long i, coord
            int tmpsum = 0

        pview [:] = 0
        for seg in segments:
            for i in range(seg.start,seg.end):
                coord = self.c_get_segmentchain_coordinate(i,False)
                pview[coord] = 1
                tmpsum += 1

        self._mask_segments = segments
        self.masked_length = self.length - tmpsum
        self._position_mask = pmask

        return True

    def get_masks(self):
        """Return masked positions as a list of |GenomicSegments|
        
        Returns
        -------
        list
            list of |GenomicSegments| representing masked positions
        
        See also
        --------
        SegmentChain.get_masks_as_segmentchain
        
        SegmentChain.add_masks
        
        SegmentChain.reset_masks
        """
        warn("SegmentChain.get_masks() is deprecated and will soon be removed from `plastid`. Use SegmentChain.mask_segments instead",UserWarning)
        if self._mask_segments is None:
            return []
        return self._mask_segments
    
    def get_masks_as_segmentchain(self):
        """Return masked positions as a |SegmentChain|
        
        Returns
        -------
        |SegmentChain|
            Masked positions
        
        See also
        --------
        SegmentChain.get_masks
        
        SegmentChain.add_masks

        SegmentChain.reset_masks
        """
        if self._mask_segments is None:
            return SegmentChain()

        return SegmentChain(*self._mask_segments)
    
    cdef void c_reset_masks(self):
        """Resets all masks and hashes that depend on them to empty values
        """
        if self._mask_segments is None:
            pass
            
        self._position_mask = None
        self._mask_segments = None 
        self.masked_length = self.length

    def reset_masks(self):
        """Removes masks added by :py:meth:`add_masks`

        See also
        --------
        SegmentChain.add_masks
        """
        self.c_reset_masks()

    def get_junctions(self):
        """Returns a list of |GenomicSegments| representing spaces
        between the |GenomicSegments| in `self` In the case of a transcript,
        these would represent introns. In the case of an alignment, these
        would represent gaps in the query compared to the reference.
        
        Returns
        -------
        list
            List of |GenomicSegments| covering spaces between the intervals in `self`
            (e.g. introns in the case of a transcript, or gaps in the case of
            an alignment)
        """
        cdef:
            GenomicSegment seg1, seg2
            list juncs = []
            GenomicSegment span = self.spanning_segment
            str strand = self.strand
            str chrom = self.chrom
            int i

        for i in range(len(self._segments)-1):
            seg1, seg2 = self._segments[i], self._segments[i+1]
            juncs.append(GenomicSegment(chrom,
                                        seg1.end,
                                        seg2.start,
                                        strand))
        return juncs
    
    # TODO: test optimality
    def as_gff3(self, str feature_type=None, bint escape=True, list excludes=None):
        """Format `self` as a line of `GFF3`_ output.
        
        Because `GFF3`_ files permit many schemas of parent-child hierarchy,
        and in order to reduce confusion and overhead, attempts to export
        a multi-interval |SegmentChain| will raise an AttributeError.
        
        Instead, users may export the individual features from which the
        multi-interval |SegmentChain| was constructed, or construct features
        for them, setting *ID*, *Parent*, and *type* attributes following
        their own conventions.

         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the `type` attribute of `self.attr`
        
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: `True`)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: `[]`)
        
        
        Returns
        -------
        str
            Line of `GFF3`_-formatted text
        
            
        Raises
        ------
        AttributeError
            if the |SegmentChain| has multiple intervals
            
        Notes
        -----
        Columns of `GFF3`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes                       
            ======== =========

        For further information, see
            - `GFF3 file format specification <http://www.sequenceontology.org/gff3.shtml>`_
            - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_
            - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_            
        """
        cdef int length = len(self._segments)

        if length == 0: # empty SegmentChain
            return ""
        elif length > 1:
            raise AttributeError("Attempted export of multi-interval %s" % self.__class__)

        cdef:
            dict gff_attr
            list always_excluded, ltmp
            GenomicSegment segment
        
        excludes = [] if excludes is None else excludes
        gff_attr = copy.deepcopy(self.attr)
        feature_type = self.attr["type"] if feature_type is None else feature_type
        
        always_excluded = ["source",
                           "score",
                           "phase",
                           "cds_genome_start",
                           "cds_genome_end",
                           "thickstart",
                           "thickend",
                           "type",
                           "_bedx_column_order"]

        for segment in self._segments:
            ltmp = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GFF3_tokens(gff_attr,
                                     excludes=always_excluded+excludes,
                                     escape=True)]

        return "\t".join(ltmp) + "\n"
    
    # TODO: optimize
    def as_gtf(self, str feature_type=None, bint escape=True, list excludes=None):
        """Format |SegmentChain| as a block of `GTF2`_ output.
        
        The `frame` or `phase` attribute (`GTF2`_ column 8) is valid only for `'CDS'`
        features, and, if not present in `self.attr`, is calculated assuming
        the |SegmentChain| contains the entire coding region. If the |SegmentChain|
        contains multiple intervals, the `frame` or `phase` attribute will
        *always* be recalculated.
        
        All attributes in `self.attr`, except those created upon import,
        will be propagated to all of the features that are generated.
        
        Parameters
        ----------
        feature_type : str
            If not None, overrides the "type" attribute of ``self.attr``
        
        escape : bool, optional
            Escape tokens in column 9 of GTF output (Default: True)
        
        excludes : list, optional
            List of attribute key names to exclude from column 8
            (Default: *[]*)
        
        Returns
        -------
        str
            Block of GTF2-formatted text
        
        
        Notes
        -----
        `gene_id` and `transcript_id` are required
            The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
            that attributes `gene_id` and `transcript_id` be defined. If these
            are not present in `self.attr`, their values will be guessed 
            following the rules in :py:meth:`SegmentChain.get_gene` and 
            :py:meth:`SegmentChain.get_name`, respectively.
        
        Beware of attribute loss
            To save memory, only the attributes shared by all of the individual
            sub-features (e.g. exons) that were used to assemble this |Transcript|
            have been stored in `self.attr`. This means that upon re-export to `GTF2`_,
            these sub-features will be lacking any attributes that were specific
            to them individually. Formally, this is compliant with the 
            `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
            explicitly that only the attributes `gene_id` and `transcript_id`
            are supported.
            
        Columns of `GTF2`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes. "gene_id" and "transcript_id" are required                        
            ======== =========
        
        For more info
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_           
        """
        if len(self) == 0:
            return ""

        cdef:
            dict gtf_attr
            dict attr = self.attr
            list ltmp1, ltmp2, always_excluded
        
        excludes = [] if excludes is None else excludes
        gtf_attr = copy.deepcopy(attr)
        gtf_attr["transcript_id"] = attr.get("transcript_id",self.get_name())
        gtf_attr["gene_id"]       = attr.get("gene_id",self.get_gene())
        feature_type = attr["type"] if feature_type is None else feature_type
        
        ltmp1 = []
        
        always_excluded = ["source",
                           "Parent",
                           "score",
                           "phase",
                           "cds_genome_start",
                           "cds_genome_end",
                           "thickstart",
                           "thickend",
                           "type",
                           "color",
                           "_bedx_column_order"]
        
        for segment in self:
            ltmp2 = self._get_8_gff_columns(segment,feature_type) +\
                   [make_GTF2_tokens(gtf_attr,
                                     excludes=always_excluded+excludes,
                                     escape=escape)]
            
            ltmp1.append("\t".join(ltmp2))            

        return "\n".join(ltmp1) + "\n"
    
    # TODO: test optimality
    def _get_8_gff_columns(self, GenomicSegment segment, str feature_type):
        """Format columns 1-8 of GFF/GTF2/GFF3 files.
        
        Parameters
        ----------
        segment : |GenomicSegment|
            Segment to export
        
        feature_type : str
            Type of feature (for column 3 of output)
        
        
        Notes
        ------
        Columns of GFF files are as follows:
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if CDS) 
                9     Attributes. Formatting depends on flavor of GFF                      
            ======== =========        
        """
        cdef:
            str chrom  = segment.chrom
            str strand = segment.strand
            str phase  = "."
            long segment_start = segment.start
            long segment_end   = segment.end
            long new_segment_start
            list ltmp
            dict attr = self.attr

        if feature_type == "CDS":
            # use phase/frame if known for length-1 features
            # called "phase" in GFF3 conventions; "frame" in GTF2
            if len(self._segments) == 1 and ("phase" in attr or "frame" in attr):
                phase = str(attr.get("phase",attr.get("frame")))
            # otherwise calculate
            else:
                new_segment_start = self.get_segmentchain_coordinate(chrom,segment_start,strand,stranded=False)
                phase = str((3 - (new_segment_start % 3)) % 3)
        
        ltmp = [chrom,
                self.attr.get("source","."),
                feature_type,
                str(segment_start + 1),
                str(segment_end),
                str(self.attr.get("score",".")),
                strand,
                phase]
        
        return ltmp
    
    # TODO: optimize
    def as_bed(self, thickstart=None, thickend=None, as_int=True, color=None, extra_columns=None, empty_value=""):
        """Format |SegmentChain| as a string of BED12[+X] output.
        
        If the |SegmentChain| was imported as a `BED`_ file with extra columns,
        these will be output in the same order, after the `BED`_ columns.

        Parameters
        ----------
        thickstart : int or `None`, optional
            If not `None`, overrides the genome coordinate that starts thick
            plotting in genome browser found in `self.attr['thickstart']`
    
        thickend : int or None, optional
            If not None, overrides the genome coordinate that stops
            thick plotting in genome browser found in `self.attr['thickend']`

        as_int : bool, optional
            Force `score` to integer (Default: `True`)
   
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in `self.attr['color']`
    
        extra_columns : None or list-like, optional
            If `None`, and the |SegmentChain| was imported using the `extra_columns`
            keyword of :meth:`~plastid.genomics.roitools.SegmentChain.from_bed`,
            the |SegmentChain| will be exported in BED 12+X format, in which
            extra columns are in the same order as they were upon import. If no extra columns
            were present, the |SegmentChain| will be exported as a BED12 line.

            If a list of attribute names, these attributes will be exported as
            extra columns in order, overriding whatever happened upon import. 
            If an attribute name is not in the `attr` dict of the |SegmentChain|,
            it will be exported with the value of `empty_value`

            If an empty list, no extra columns will be exported; the |SegmentChain|
            will be formatted as a BED12 line.
        
        empty_value : str, optional
            Value to export for `extra_columns` that are not defined (Default: "")


        Returns
        -------
        str 
            Line of BED12[+X]-formatted text


        Notes
        -----
        BED12 columns are as follows:
            ======== =========
            Column   Contains
            ======== =========
               1     Contig or chromosome
               2     Start of first block in feature (0-indexed)
               3     End of last block in feature (half-open)
               4     Feature name
               5     Feature score
               6     Strand
               7     thickstart (in chromosomal coordinates)
               8     thickend (in chromosomal coordinates)
               9     Feature color as RGB tuple
               10    Number of blocks in feature
               11    Block lengths
               12    Block starts, relative to start of first block
            ======== =========

        For more details
            See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        cdef:
            list ltmp
            GenomicSegment span = self.spanning_segment
            long spanstart = span.start
            
        if len(self) > 0:
            score = self.attr.get("score",0)
            try:
                score = float(score)
                if as_int is True:
                    score = int(round(score))
            except ValueError:
                score = 0
            except TypeError:
                score = 0
            
            try:
                color = self.attr.get("color","#000000")
                if isinstance(color,str):
                    color = "%s,%s,%s" % tuple(get_rgb255(color))
                else:
                    # is a tuple or sequence
                    if all(isinstance(X,float) for X in color) and \
                       all(X < 1 for X in color) and \
                       all(X > 0 for X in color):
                        color = "%s,%s,%s" % (int(round(255*X)) for X in color)
                    else:
                        color = "%s,%s,%s" % tuple(color)

            except ValueError:
                color = self.attr.get("color","0,0,0") if color is None else color
            
            thickstart = self.attr.get("thickstart",spanstart) if thickstart is None else thickstart
            thickend   = self.attr.get("thickend",spanstart)   if thickend   is None else thickend
            
            ltmp = [span.chrom,
                    spanstart,
                    span.end,
                    self.get_name(),
                    score,
                    span.strand,
                    thickstart,
                    thickend,
                    color,
                    len(self),
                    ",".join([str(len(X)) for X in self]) + ",",
                    ",".join([str(X.start - self[0].start) for X in self]) + ","
                   ]            

            if extra_columns is None:
                extra_columns = self.attr.get("_bedx_column_order",[])

            if len(extra_columns) > 0:
                ltmp.extend([self.attr.get(X,empty_value) for X in extra_columns])
            
            return "\t".join([str(X) for X in ltmp]) + "\n"
        else:
            # SegmentChain with no intervals
            return ""
    
    def as_psl(self):
        """Formats |SegmentChain| as `PSL`_ (blat) output.
        
        Notes
        -----
        This will raise an :py:class:`AttributeError` unless the following
        keys are present and defined in `self.attr`, corresponding to the
        columns of a `PSL`_ file:
        
            ======  ===================================
            Column  Key
            ======  ===================================
                1   ``match_length``
                2   ``mismatches``
                3   ``rep_matches``
                4   ``N``
                5   ``query_gap_count``
                6   ``query_gap_bases``
                7   ``target_gap_count``
                8   ``target_gap_bases``
                9   ``strand``
                10  ``query_name``
                11  ``query_length``
                12  ``query_start``
                13  ``query_end``
                14  ``target_name``
                15  ``target_length``
                16  ``target_start``
                17  ``target_end``
                19  ``q_starts`` : list of integers
                20  ``l_starts`` : list of integers
            ======  ===================================
        
        These keys are defined only if the |SegmentChain| was created by
        :py:meth:`SegmentChain.from_psl`, or if the user has defined them.
        
        See the `PSL spec <http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html>`_
        for more information.
        
        
        Returns
        -------
        str
            PSL-representation of BLAT alignment

        
        Raises
        ------
        AttributeError
            If not all of the attributes listed above are defined
        """
        cdef:
            dict attr = self.attr
            list ltmp
            str block_sizes, q_starts, t_starts
        try:
            ltmp = [
                attr["match_length"],
                attr["mismatches"],
                attr["rep_matches"],
                attr["N"],
                attr["query_gap_count"],
                attr["query_gap_bases"],
                attr["target_gap_count"],
                attr["target_gap_bases"],
                attr["strand"],
                attr["query_name"],
                attr["query_length"],
                attr["query_start"],
                attr["query_end"],
                attr["target_name"],
                attr["target_length"],
                attr["target_start"],
                attr["target_end"],
                len(self),
            ]
   
            block_sizes = ",".join([str(len(X)) for X in self]) + ","
            q_starts = ",".join([str(X) for X in self.attr["q_starts"]]) + ","
            t_starts = ",".join([str(X) for X in self.attr["t_starts"]]) + ","
            
            ltmp.append(block_sizes)
            ltmp.append(q_starts)
            ltmp.append(t_starts)
            return "\t".join(str(X) for X in ltmp) + "\n"
        except KeyError:
            raise AttributeError("SegmentChains only support PSL output if all PSL attributes are defined in self.attr: match_length, mismatches, rep_matches, N, query_gap_count, query_gap_bases, strand, query_length, query_start, query_end, target_name, target_length, target_start, target_end")
    
    def get_segmentchain_coordinate(self, str chrom, long genomic_x, str strand, bint stranded = True):
        """Finds the |SegmentChain| coordinate corresponding to a genomic position
        
        Parameters
        ----------
        chrom : str
            Chromosome name
            
        genomic_x : int
            coordinate, in genomic space
            
        strand : str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
            
        stranded : bool, optional
            If `True`, coordinates are given in stranded space
            (i.e. from 5' end of chain, as one might expect for a transcript).
            If `False`, coordinates are given from the left end of `self`,
            regardless of strand. (Default: `True`)
        
        
        Returns
        -------
        int
            Position in |SegmentChain|
            
        Raises
        ------
        KeyError
            if position outside bounds of |SegmentChain|
        """
        cdef long retval
        if chrom != self.chrom:
            raise ValueError("get_segmentchain_coordinate: query chromosome '%s' does not match chain '%s'" % chrom, self)
        if strand != self.strand:
            raise ValueError("get_segmentchain_coordinate: query strand '%s' does not match chain '%s'" % strand, self)

        retval = self.c_get_segmentchain_coordinate(genomic_x,stranded)
        if retval == -1:
            raise KeyError("SegmentChain.get_segmentchain_coordinate: Genomic position %s not in SegmentChain %s." % (genomic_x,self))

        return retval

    cdef long c_get_segmentchain_coordinate(self, long genomic_x, bint stranded) except -1:
        """Finds the |SegmentChain| coordinate corresponding to a genomic position
        
        Parameters
        ----------
        genomic_x : int
            coordinate, in genomic space
            
        stranded : bool, optional
            If `True`, coordinates are given in stranded space
            (i.e. from 5' end of chain, as one might expect for a transcript).
            If `False`, coordinates are given from the left end of `self`,
            regardless of strand. (Default: `True`)
        
        
        Returns
        -------
        int
            Position in |SegmentChain|
            
        Raises
        ------
        KeyError
            if position outside bounds of |SegmentChain|
        """
        cdef:
            int i = 0
            int num_segs = len(self._segments)
            long cumlength = 0
            long retval
            long chrom_end
            GenomicSegment span = self.spanning_segment
            GenomicSegment seg
            str msg

        msg = "SegmentChain.get_segmentchain_coordinate: genomic position '%s' is not in chain '%s'.\n" % (genomic_x,self)

        if genomic_x < span.start:
            raise KeyError(msg)

        while i < num_segs:
            seg = self._segments[i]
            chrom_end = seg.end
            cumlength += len(seg)
            if genomic_x < chrom_end:
                if genomic_x >= self._segments[i].start:
                    retval = cumlength  - chrom_end + genomic_x
                    if span.c_strand == reverse_strand and stranded == True:
                        retval = self.length - retval - 1
                    return retval
                # because segments are sorted, if < end but not >= start,
                # key must be outside bounds of chain
                raise KeyError(msg) 
            i += 1

        raise KeyError(msg)

    def get_genomic_coordinate(self,x,stranded=True):
        """Finds genomic coordinate corresponding to position `x` in `self`
        
        Parameters
        ----------
        x : int
            position of interest, relative to |SegmentChain|
            
        stranded : bool, optional
            If `True`, `x` is assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            coordinates assumed to be counted the left end of the `self`,
            regardless of the strand of `self`. (Default: `True`)
        
                             
        Returns
        -------
        str 
            Chromosome name
        
        long
            Genomic cordinate corresponding to position `x`
        
        str
            Chromosome strand (`'+'`, `'-'`, or `'.'`)
        
        
        Raises
        ------
        IndexError
            if `x` is outside the bounds of the |SegmentChain|
        """
        cdef long new_x = self.c_get_genomic_coordinate(x, stranded)
        if new_x == -1:
            raise RuntimeError("Cannot fetch coordinate %s from SegmentChain %s (length %s)" % (x,self,self.length))

        return self.chrom, new_x, self.strand

    cdef long c_get_genomic_coordinate(self, long x, bint stranded) except -1:
        """Finds genomic coordinate corresponding to position `x` in `self`
        
        Parameters
        ----------
        x : int
            position of interest, relative to |SegmentChain|
            
        stranded : bool
            If `True`, `x` is assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            coordinates assumed to be counted the left end of the `self`,
            regardless of the strand of `self`.

        .. note::

           Unlike :meth:`SegmentChain.get_genomic_coordinate`, this only returns
           a `long`, as opposed to a tuple of (chrom_name, position, chrom_strand)
           

        Returns
        -------
        long
            Genomic coordinate corresponding to position `x`


        See also
        --------
        SegmentChain.get_genomic_coordinate
        """
        cdef:
            int i = 0
            int num_segs = len(self._segments)
            long retval
            long lastlength = 0
            long cumlength = 0
            long length = self.length
            GenomicSegment seg
            str msg

        msg = "Position %s is outside bounds [0,%s) of SegmentChain '%s'" % \
              (x, length, self.get_name())
        if x < 0 or x >= length:
            raise IndexError(msg)

        if self.c_strand == reverse_strand and stranded == True:
            x = length - x - 1

        if self._position_hash is None:
            while i < num_segs:
                seg = self._segments[i]
                cumlength += len(seg)
                if x < cumlength:
                    if x >= lastlength:
                        retval = x - lastlength + seg.start
                        
                        return retval
                    raise IndexError(msg)

                lastlength = cumlength
                i += 1
        else:
            return self._position_hash[x]
        
        raise IndexError(msg)

    def get_subchain(self, long start, long end, bint stranded=True, **extra_attr):
        """Retrieves a sub-|SegmentChain| corresponding a range of positions
        specified in coordinates relative this |SegmentChain|. Attributes in
        `self.attr` are copied to the child SegmentChain, with the exception
        of `ID`, to which the suffix `'subchain'` is appended.
        
        Parameters
        ----------
        start : int
            position of interest in SegmentChain coordinates, 0-indexed
            
        end : int
            position of interest in SegmentChain coordinates, 0-indexed 
            and half-open
            
        stranded : bool, optional
            If `True`, `start` and `end` are assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            they assumed to be counted the left end of the `self`,
            regardless of the strand of `self`. (Default: `True`)

        extra_attr : keyword arguments
            Values that will be included in the subchain's `attr` dict.
            These can be used to overwrite values already present.
                          
        Returns
        -------
        |SegmentChain|
            covering parent chain positions `start` to `end` of `self`
        
        
        Raises
        ------
        IndexError
            if `start` or `end` is outside the bounds of the |SegmentChain|

        TypeError
            if `start` or `end` is None
        """
        cdef:
            SegmentChain chain = self.c_get_subchain(start,end,stranded)
            dict old_attr = copy.deepcopy(self.attr)

        old_attr.update(extra_attr)
        chain.attr = old_attr
        chain.attr["ID"] = "%s_subchain" % self.get_name()
        return chain

    cdef SegmentChain c_get_subchain(self, long start, long end, bint stranded):
        """Similar to :meth:`SegmentChain.get_subchain` but does not copy `attr` dict
        
        Parameters
        ----------
        start : long
            position of interest in SegmentChain coordinates, 0-indexed
            
        end : long
            position of interest in SegmentChain coordinates, 0-indexed 
            and half-open
            
        stranded : bool
            If `True`, `start` and `end` are assumed to be in stranded space (i.e. counted from
            5' end of chain, as one might expect for a transcript). If `False`,
            they assumed to be counted the left end of the `self`,
            regardless of the strand of `self`.

        Returns
        -------
        |SegmentChain|
            covering parent chain positions `start` to `end`

        See also
        --------
        SegmentChain.get_subchain
        """
        cdef:
            SegmentChain chain = SegmentChain()
            long length = self.length
            long tmp
            list segs

        if start == end: # this is a special case which we need to account for
            return SegmentChain()

        if start is None:
            raise TypeError('start coordinate supplied is None. Expected int')
        elif end is None:
            raise TypeError('end coordinate supplied is None. Expected int')

        if stranded == True and self.c_strand == reverse_strand:
            tmp = end
            end   = length - start 
            start = length - tmp 
            
        positions = self._get_position_hash()[start:end].tolist()
        segs = positionlist_to_segments(self.chrom,self.strand,positions)
        chain._set_segments(segs)

        return chain

    def get_counts(self,ga,stranded=True):
        """Return list of counts or values drawn from `ga` at each position in `self`
        
        Parameters
        ----------
        ga : GenomeArray from which to fetch counts
            
        stranded : bool, optional
            If `True` and `self` is on the minus strand,
            count order will be reversed relative to genome so that the
            array positions march from the 5' to 3' end of the chain.
            (Default: `True`)
            
            
        Returns
        -------
        numpy.ndarray
            Array of counts from `ga` covering `self`
        """
        cdef:
            GenomicSegment  seg
            list segments = self._segments
            numpy.ndarray   count_array
            list count_arrays
            list dims
            int i, j, n

        if len(self) == 0:
            warn("%s is a zero-length SegmentChain. Returning 0-length count vector." % self.get_name(),DataWarning)
            return numpy.array([],dtype=float)

#        count_array = numpy.concatenate([ga.get(X,roi_order=False) for X in self._segments],axis=-1)
        
        # code block below is ugly compared to numpy.concatenate, but this is actually
        # faster. must be something weird in implementation
        count_arrays = [ga.get(X,roi_order=False) for X in segments]
        dims         = list(count_arrays[0].shape)
        dims[-1]     = self.length
        count_array  = numpy.empty(dims,dtype=float)
         
        i = 0
        for n,seg in enumerate(segments):
            j = i + len(seg)
            count_array[...,i:j] = count_arrays[n]
            i = j

        if self.c_strand == reverse_strand and stranded is True:
            count_array = count_array[...,::-1]
            
        return count_array

    def get_masked_counts(self,ga,stranded=True,copy=False):
        """Return counts covering `self` in dataset `gnd` as a masked array, in transcript 
        coordinates. Positions masked by :py:meth:`SegmentChain.add_mask` 
        will be masked in the array
        
        Parameters
        ----------
        gnd : non-abstract subclass of |AbstractGenomeArray|
            GenomeArray from which to fetch counts
            
        stranded : bool, optional
            If true and the |SegmentChain| is on the minus strand,
            count order will be reversed relative to genome so that the
            array positions march from the 5' to 3' end of the chain.
            (Default: `True`)

        copy : bool, optional
            If `False` (default) returns a view of the data; so changing
            values in the view changes the values in the |GenomeArray|
            if it is mutable. If `True`, a copy is returned instead.
            
            
        Returns
        -------
        :py:class:`numpy.ma.masked_array`
        """
        cdef:
            numpy.ndarray counts = self.get_counts(ga)
            numpy.ndarray mask, m

        if self._position_mask is None:
            mask = numpy.zeros_like(counts)
        else:
            # need to use intc bc self._position_mask is an array.array
            m = numpy.frombuffer(self._position_mask,dtype=numpy.intc)
            if self.c_strand == reverse_strand:
                m = m[::-1]
            
            mask = numpy.empty_like(counts)
            mask[...,:] = m

        return MaskedArray(counts,mask=mask.astype(bool),copy=copy)
        
    def get_sequence(self,genome,stranded=True):
        """Return spliced genomic sequence of |SegmentChain| as a string
        
        Parameters
        ----------
        genome : dict or :class:`twobitreader.TwoBitFile`
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If `True` and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (Default: True)
            
            
        Returns
        -------
        str
            Nucleotide sequence of the |SegmentChain| extracted from `genome`
        """
        cdef:
            list ltmp
            str stmp

        if len(self._segments) == 0:
            warn("%s is a zero-length SegmentChain. Returning empty sequence." % self.get_name(),DataWarning)
            return ""

        else:
            chromseq = genome[self.spanning_segment.chrom]
            ltmp = [chromseq[X.start:X.end] for X in self]
            stmp = "".join([str(X.seq) if isinstance(X,SeqRecord) else X for X in ltmp])

            if self.strand == "-"  and stranded == True:
                stmp = str(Seq(stmp,generic_dna).reverse_complement())
            
        return stmp
    
    def get_fasta(self,genome,stranded=True):
        """Formats sequence of SegmentChain as FASTA output
        
        Parameters
        ----------
        genome : dict or :class:`twobitreader.TwoBitFile`
            Dictionary mapping chromosome names to sequences.
            Sequences may be strings, string-like, or :py:class:`Bio.Seq.SeqRecord` objects
       
        stranded : bool
            If `True` and the |SegmentChain| is on the minus strand,
            sequence will be reverse-complemented (Default: True)

            
        Returns
        -------
        str
            FASTA-formatted seuqence of |SegmentChain| extracted from `genome`
        """
        return ">%s\n%s\n" % (self.get_name(),self.get_sequence(genome,stranded=stranded))

    @staticmethod
    def from_str(str inp):
        """Create a |SegmentChain| from a string formatted by :py:meth:`SegmentChain.__str__`:
           
            `chrom:start-end^start-end(strand)`
           
        where '^' indicates a splice junction between regions specified
        by `start` and `end` and `strand` is '+', '-', or '.'. Coordinates are
        0-indexed and half-open.


        Parameters
        ----------
        inp : str
			String formatted in manner of :py:meth:`SegmentChain.__str__`
          
          
        Returns
        -------
        |SegmentChain|
        """
        cdef:
            str chrom, middle, strand, sstart, send, piece
            list segs
            long start, end

        if inp in ("na","nan","None:(None)","None","none",None) or isinstance(inp,float) and numpy.isnan(inp):
            return SegmentChain()
        else:
            chrom,middle,strand = ivcpat.search(inp).groups()
            segs = []
            for piece in middle.split("^"):
                sstart,send = piece.split("-")
                start = <long>long(sstart)
                end = <long>long(send)
                segs.append(GenomicSegment(chrom,start,end,strand))
            return SegmentChain(*segs)
        
    # at each call to from_bed (e.g. from a BigBedReader)
    @staticmethod
    def from_bed(str line, object extra_columns=0):
        """Create a |SegmentChain| from a line from a `BED`_ file.
        The `BED`_ line may contain 4 to 12 columns, per the specification.
        These will be auto-detected and parsed appropriately.
        
        See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        for more details.

        Parameters
        ----------
        line
            Line from a `BED`_ file, containing 4 or more columns

        extra_columns: int or list optional
            Extra, non-BED columns in :term:`Extended BED`_ format file corresponding
            to feature attributes. This is common in `ENCODE`_-specific `BED`_ variants.
            
            if `extra-columns` is:
            
              - an :class:`int`: it is taken to be the
                number of attribute columns. Attributes will be stored in
                the `attr` dictionary of the |SegmentChain|, under names like
                `custom0`, `custom1`, ... , `customN`.

              - a :class:`list` of :class:`str`, it is taken to be the names
                of the attribute columns, in order, from left to right in the file.
                In this case, attributes in extra columns will be stored under
                their respective names in the `attr` dict.

              - a :class:`list` of :class:`tuple`, each tuple is taken
                to be a pair of `(attribute_name, formatter_func)`. In this case,
                the value of `attribute_name` in the `attr` dict of the |SegmentChain|
                will be set to `formatter_func(column_value)`.
            
            (Default: 0)

        Returns
        -------
        |SegmentChain|
        """
        cdef:
            dict         attr     = get_attr_from_bed(line,extra_columns=extra_columns)
            list         segments = get_segments_from_bed(attr)
            SegmentChain chain    = SegmentChain()
            
        chain._set_segments(segments)
        chain.attr.update(attr)
        return chain

        
    @staticmethod
    def from_psl(psl_line):
        """Create a |SegmentChain| from a line from a `PSL`_ (BLAT) file

        See the `PSL spec <http://pombe.nci.nih.gov/genome/goldenPath/help/blatSpec.html>`_
        
        Parameters
        ----------
        psl_line : str
            Line from a `PSL`_ file

        Returns
        -------
        |SegmentChain|
        """
        cdef:
            list items = psl_line.strip().split("\t")        
            list segs = []
            dict attr = {}
            long t_start, block_size
            str st, sb
            GenomicSegment seg
            list block_sizes, sq_starts, st_starts
            list t_starts = []
            list q_starts = []

        attr["type"]             = "alignment"
        attr["query_name"]       = items[9]
        attr["match_length"]     = int(items[0])
        attr["mismatches"]       = int(items[1])
        attr["rep_matches"]      = int(items[2])
        attr["N"]                = int(items[3])
        attr["query_gap_count"]  = int(items[4])
        attr["query_gap_bases"]  = int(items[5])
        attr["target_gap_count"] = int(items[6])
        attr["target_gap_bases"] = int(items[7])
        attr["strand"]           = items[8]
        attr["query_length"]     = int(items[10])
        attr["query_start"]      = int(items[11])
        attr["query_end"]        = int(items[12])
        attr["target_name"]      = items[13]
        attr["target_length"]    = int(items[14])
        attr["target_start"]     = int(items[15])
        attr["target_end"]       = int(items[16])
        attr["ID"]               = attr["query_name"]
        #block_count           = int(items[17])

        block_sizes = items[18].strip(",").split(",")
        sq_starts    = items[19].strip(",").split(",")
        st_starts    = items[20].strip(",").split(",")        

        for st, sb, sq in zip(st_starts,block_sizes,sq_starts):
            t_start = long(st)
            block_size = long(sb)
            t_starts.append(t_start)
            q_starts.append(long(sq))
            seg = GenomicSegment(attr["target_name"],
                                 t_start,
                                 t_start + block_size,
                                 attr["strand"])
            segs.append(seg)
        
        attr["t_starts"] = t_starts 
        attr["q_starts"] = q_starts

        return SegmentChain(*segs,**attr)        


 
cdef class Transcript(SegmentChain):
    """
    Transcript(*segments,**attributes)
    
    Subclass of |SegmentChain| specifically for RNA transcripts.
    In addition to coordinate-conversion, count fetching, sequence fetching,
    and various other methods inherited from |SegmentChain|, |Transcript|
    provides convenience methods for fetching sub-chains corresponding to 
    CDS features, 5' UTRs, and 3' UTRs.


    Parameters
    ----------
    *segments : |GenomicSegment|
        0 or more GenomicSegments on the same strand
    
    **attr : keyword arguments
        Arbitrary attributes, including, for example:
    
        ====================    ============================================
        **Attribute**           **Description**
        ====================    ============================================
        ``cds_genome_start``    Location of CDS start, in genomic coordinates

        ``cds_genome_start``    Location of CDS end, in genomic coordinates
        
        ``ID``                  A unique ID for the |SegmentChain|.

        ``transcript_id``       A transcript ID used for `GTF2`_ export

        ``gene_id``             A gene ID used for `GTF2`_ export
        ====================    ============================================    




    Attributes
    ----------
    cds_genome_start : int or None
        Leftmost position in genomic coordinates of coding region, 0-indexed

    cds_genome_end : int or None
        Rightmost position in genomic coordinates of coding region, 0-indexed
        and half-open

    cds_start : int or None
        Stranded position relative to 5' end of transcript at which coding region starts
        (note: for minus-strand features this will be higher in genomic
        coordinates than `cds_end`).

    cds_end : int or None
        Stranded position relative to 5' end of transcript at which coding region ends
        (note: for minus-strand features this will be lower in genomic coordinates
        than `cds_start`).

    spanning_segment : |GenomicSegment|
        A GenomicSegment spanning the endpoints of the Transcript

    strand : str
        The chromosome strand ('+', '-', or '.')

    chrom : str
        The chromosome name
        
    segments : list
        List of |GenomicSegments| in the chain
    
    mask_segments : list
        List of |GenomicSegments| that mask positions in the chain.

    attr : dict
        Miscellaneous attributes
    """
    
    def __cinit__(self,*segments,type=None,**attr):
        """Create a |Transcript|
        
        Parameters
        ----------
        *segments : |GenomicSegment|
            0 or more |GenomicSegments| (exons)

        **attr : dict
            keyword attributes

        attr["cds_genome_start"] : int or None
            genome coordinate of CDS start, if any
                         
        attr["cds_genome_end"] : int or None
            genome coordinate of CDS end, if any
    
        attr["type"] : str
            If provided, a feature type used for GTF2/GFF3 export
            Otherwise, set to "mRNA"
        
        attr["ID"] : str
            If provided, a unique ID for the |Transcript|.
            Otherwise, generated from genomic coordinates
        
        attr["transcript_id"] : str
            If provided, a transcript_id used for `GTF2`_ export.
            Otherwise, generated from genomic coordinates.
        
        attr["gene_id"] : str
            If provided, a gene_id used for `GTF2`_ export
            Otherwise, generated from genomic coordinates.
        """
        if type == None:
            self.attr["type"] = "mRNA"
        
        gstart = attr.get("cds_genome_start",None) 
        gend   = attr.get("cds_genome_end"  ,None)
        self.cds_genome_start = gstart
        self.cds_genome_end   = gend

        if gstart is not None and gend is not None:
            self._update_cds()
        else:
            self.cds_genome_start = None
            self.cds_genome_end = None
            self.cds_start = None
            self.cds_end = None
 
    def __copy__(self):  # copy info and segments; shallow copy attr
        chain2 = Transcript()
        chain2._set_segments(self._segments)
        if self._mask_segments is not None:
            chain2._set_masks(self._mask_segments)
        chain2.cds_genome_start = self.cds_genome_start
        chain2.cds_genome_end = self.cds_genome_end
        if chain2.cds_genome_start is not None:
            chain2._update_from_cds_genome_start()
        if chain2.cds_genome_end is not None:
            chain2._update_from_cds_genome_end()
        chain2.attr = self.attr
        return chain2

    def __deepcopy__(self,memo): # deep copy everything
        chain2 = Transcript()
        chain2._set_segments(copy.deepcopy(self._segments))
        if self._mask_segments is not None:
            chain2._set_masks(copy.deepcopy(self._mask_segments))
        chain2.cds_genome_start = self.cds_genome_start
        chain2.cds_genome_end = self.cds_genome_end
        if chain2.cds_genome_start is not None:
            chain2._update_from_cds_genome_start()
        if chain2.cds_genome_end is not None:
            chain2._update_from_cds_genome_end()

        chain2.attr = copy.deepcopy(self.attr)
        return chain2

    def __reduce__(self):
        return (Transcript,tuple(),self.__getstate__())

    def __getstate__(self): # pickle state
        cdef:
            list segstrs  = [str(X) for X in self._segments]
            list maskstrs = []
            dict attr = self.attr
        if self._mask_segments is not None:
            maskstrs = [str(X) for X in self._mask_segments]

        return (segstrs, maskstrs, attr, self.cds_genome_start, self.cds_genome_end)

    def __setstate__(self,state): # revive state from pickling
        cdef:
            list segs  = [GenomicSegment.from_str(X) for X in state[0]]
            list masks = [GenomicSegment.from_str(X) for X in state[1]]
            dict attr = state[2]
            object gstart = state[3]
            object gend = state[4]

        self.attr = attr
        self._set_segments(segs)
        if len(masks) > 0:
            self._set_masks(masks)
        self.cds_genome_start = gstart
        self.cds_genome_end = gend
        if gstart is not None and gend is not None:
            self._update_cds()
        else:
            self.cds_start = None
            self.cds_end = None

    def __richcmp__(self, object other, int cmpval):
        """Test whether `self` and `other` are equal or unequal. Equality is defined as
        identity of positions, chromosomes, and strands. Two |SegmentChain| with
        zero intervals, by convention, are not equal.
           
        Parameters
        ----------
        other : |SegmentChain| or |GenomicSegment|
        	Query feature
        
        Returns
        -------
        bool
        """
        if isinstance(other,Transcript):
            return transcript_richcmp(self,other,cmpval) == true
        elif isinstance(other,SegmentChain):
            return chain_richcmp(self,other,cmpval) == true
        elif isinstance(other,GenomicSegment):
            return chain_richcmp(self,SegmentChain(other),cmpval) == true 
        else:
            raise TypeError("SegmentChain eq/ineq/et c is only defined for other SegmentChains or GenomicSegments.")

    property cds_start:
        """Start of coding region relative to 5' end of transcript, in direction of transcript.
        Setting to None also sets `self.cds_end`, `self.cds_genome_start` and
        `self.cds_genome_end` to None
        """
        def __get__(self):
            return self.cds_start
        def __set__(self,val):
            cdef object end = self.cds_end
            if val is not None and end is not None:
                if val > end:
                    raise ValueError("Transcript '%s': cds_start (%s) must be <= cds_end (%s)" % (self,val,end))
                if val < 0:
                    raise ValueError("Transcript '%s': cds_start (%s) must be >= 0" % (self,val))
            self.cds_start = val
            if val is None:
                self.cds_end = None
                self.cds_genome_start = None
                self.cds_genome_end = None
            else:
                self._update_from_cds_start()

    property cds_end:
        """End of coding region relative to 5' end of transcript, in direction of transcript.
        Setting to None also sets `self.cds_start`, `self.cds_genome_start` and
        `self.cds_genome_end` to None
        """
        def __get__(self):
            return self.cds_end
        def __set__(self,val):
            cdef object start = self.cds_start
            if val is not None and start is not None:
                if val < start:
                    raise ValueError("Transcript '%s': cds_end (%s) must be >= cds_start (%s)" % (self,val,start))
                if val > self.length:
                    raise ValueError("Transcript '%s': cds_end (%s) must be <= self.length (%s)" % (self,val,self.length))
            self.cds_end = val
            if val is None:
                self.cds_start = None
                self.cds_genome_start = None
                self.cds_genome_end = None
            else:
                self._update_from_cds_end()

    property cds_genome_start:
        """Starting coordinate of coding region, relative to genome (i.e. leftmost;
        is start codon for forward-strand features, stop codon for reverse-strand
        features). Setting to None also sets `self.cds_start`, `self.cds_end`, and
        `self.cds_genome_end` to None
        """
        def __get__(self):
            return self.cds_genome_start
        def __set__(self, val):
            cdef object end = self.cds_genome_end
            if val is not None and end is not None:
                if val > end:
                    raise ValueError("Transcript '%s': cds_genome_start (%s) must be <= cds_genome_end (%s)" % (self,val,end))
            self.cds_genome_start = val
            if val is None:
                self.cds_genome_end = None
                self.cds_start = None
                self.cds_end = None
            else:
                self._update_from_cds_genome_start()

    property cds_genome_end:
        """Ending coordinate of coding region, relative to genome (i.e. leftmost;
        is stop codon for forward-strand features, start codon for reverse-strand
        features. Setting to None also sets `self.cds_start`, `self.cds_end`,
        and `self.cds_genome_start` to None
        """
        def __get__(self):
            return self.cds_genome_end
        def __set__(self, val):
            cdef object start = self.cds_genome_start
            if val is not None and start is not None:
                if val < start:
                    raise ValueError("Transcript '%s': cds_genome_end (%s) must be >= cds_genome_start (%s)" % (self,val,start))
            self.cds_genome_end = val
            if val is None:
                self.cds_genome_start = None
                self.cds_start = None
                self.cds_end = None
            else:
                self._update_from_cds_genome_end()

    # TODO: create explicit unit test
    cdef bint _update_cds(self) except False:
        """Generate `self.cds_start` and `self.cds_end` from `self.cds_genome_start` and `self.cds_genome_end`
        AFTER verifying that these values are not None
        """
        cdef:
            GenomicSegment span = self.spanning_segment
            Strand strand = span.c_strand
            str chrom = span.chrom
            long cds_genome_start = self.cds_genome_start
            long cds_genome_end   = self.cds_genome_end

        if strand == forward_strand:
            self.cds_start = self.c_get_segmentchain_coordinate(cds_genome_start,True)
            
            # this is in a try-catch because if the half-open cds_end coincides
            # with the end of an exon, it will not be in the end-inclusive position
            try:
                self.cds_end = self.c_get_segmentchain_coordinate(cds_genome_end,True)
            except KeyError:
                # minus one, plus one corrections because end-exclusive genome
                # position will not be in position hash if it coincides with
                # the end of any exon
                self.cds_end   = 1 + self.c_get_segmentchain_coordinate(cds_genome_end - 1,True)
        else:
            # likewise for minus-strand
            # both this adjustment and the one above for plus-strand features
            # have been thoroughly tested by examining BED files exported
            # for this purpose
            self.cds_start = self.c_get_segmentchain_coordinate(cds_genome_end - 1,  True)
            self.cds_end   = 1 + self.c_get_segmentchain_coordinate(cds_genome_start,True)

        return True

    cdef bint _update_from_cds_start(self) except False:
        cdef:
            long cds_start = <long>self.cds_start
            long [:] phash = self._get_position_hash()

        if self.spanning_segment.c_strand == forward_strand:
            self.cds_genome_start = phash[cds_start]
        else:
            self.cds_genome_end = phash[self.length - cds_start - 1] + 1 # CHECKME

        return True
    
    cdef bint _update_from_cds_end(self) except False:
        cdef:
            long cds_end = <long>self.cds_end
            long [:] phash = self._get_position_hash()

        if self.spanning_segment.c_strand == forward_strand:
            self.cds_genome_end = phash[cds_end - 1] + 1
        else:
            self.cds_genome_start = phash[self.length - cds_end]

        return True
    
    cdef bint _update_from_cds_genome_start(self) except False:
        """Generate `self.cds_start` and `self.cds_end` from `self.cds_genome_start` and `self.cds_genome_end`
        AFTER verifying that these values are not None
        """
        cdef:
            long cds_genome_start = <long>self.cds_genome_start

        if self.spanning_segment.c_strand == forward_strand:
            self.cds_start = self.c_get_segmentchain_coordinate(cds_genome_start,True)
        else:
            self.cds_end   = 1 + self.c_get_segmentchain_coordinate(cds_genome_start,True)

        return True

    cdef bint _update_from_cds_genome_end(self) except False:
        """Generate `self.cds_start` and `self.cds_end` from `self.cds_genome_start` and `self.cds_genome_end`
        AFTER verifying that these values are not None
        """
        cdef:
            GenomicSegment span = self.spanning_segment
            Strand strand = span.c_strand
            str chrom = span.chrom
            long cds_genome_start = self.cds_genome_start
            long cds_genome_end   = self.cds_genome_end

        if strand == forward_strand:
            # this is in a try-catch because if the half-open cds_end coincides
            # with the end of an exon, it will not be in the end-inclusive position
            try:
                self.cds_end = self.c_get_segmentchain_coordinate(cds_genome_end,True)
            except KeyError:
                # minus one, plus one corrections because end-exclusive genome
                # position will not be in position hash if it coincides with
                # the end of any exon
                self.cds_end   = 1 + self.c_get_segmentchain_coordinate(cds_genome_end - 1,True)
        else:
            # likewise for minus-strand
            # both this adjustment and the one above for plus-strand features
            # have been thoroughly tested by examining BED files exported
            # for this purpose
            self.cds_start = self.c_get_segmentchain_coordinate(cds_genome_end - 1,  True)

        return True

    def get_name(self):
        """Return the name of `self`, first searching through
        `self.attr` for the keys `transcript_id`, `ID`, `Name`, and `name`.
        If no value is found, :meth:`Transcript.__str__` is used.
        
        Returns
        -------
        str
            Returns in order of preference, `transcript_id`, `ID`, `Name`,
            or `name` from `self.attr`. If not found, returns ``str(self)``
        """
        cdef str name
        name = self.attr.get("transcript_id",
               self.attr.get("ID",
               self.attr.get("Name",
               self.attr.get("name",
                             str(self)))))
        return name
   
    def get_cds(self,**extra_attr):
        """Retrieve |SegmentChain| covering the coding region of `self`, including the stop codon.
        If no coding region is present, returns an empty |SegmentChain|.
        
        The following attributes are passed from `self.attr` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_CDS % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the CDS subchain's `attr` dict.
            These can be used to overwrite values already present.
        
        Returns
        -------
        |SegmentChain|
            CDS region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef:
            SegmentChain chain
            Transcript   transcript = Transcript()
            dict         attr       = {}

        if self.cds_genome_start is not None and self.cds_genome_end is not None:
            chain = self.c_get_subchain(self.cds_start,
                                        self.cds_end,
                                        True)

            transcript._set_segments(chain._segments)
            attr["gene_id"] = self.get_gene()
            attr["transcript_id"] = self.get_name()
            attr["ID"] = "%s_CDS" % self.get_name()

            transcript.attr = attr
            transcript.attr.update(extra_attr)
            transcript.cds_genome_start = self.cds_genome_start
            transcript.cds_genome_end   = self.cds_genome_end
            transcript._update_cds()
            
        return transcript
    
    def get_utr5(self,**extra_attr):
        """Retrieve sub-|SegmentChain| covering 5'UTR of `self`.
        If no coding region, returns an empty |SegmentChain|

        The following attributes are passed from `self.attr` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_5UTR" % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the 5'UTR subchain's `attr` dict.
            These can be used to overwrite values already present.


        Returns
        -------
        |SegmentChain|
            5' UTR region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef:
            SegmentChain my_segmentchain
            dict attr = {}

        if self.cds_genome_start is not None and self.cds_genome_end is not None:

            my_segmentchain = self.c_get_subchain(0,self.cds_start,True)
            attr["type"] = "5UTR"
            attr["gene_id"] = self.get_gene()
            attr["transcript_id"] = self.get_gene()
            attr["ID"] = "%s_5UTR" % self.get_name()
            my_segmentchain.attr = attr
            my_segmentchain.attr.update(extra_attr)

            return my_segmentchain
        else:
            return SegmentChain()
    
    def get_utr3(self,**extra_attr):
        """Retrieve sub-|SegmentChain| covering 3'UTR of `self`, excluding
        the stop codon. If no coding region, returns an empty |SegmentChain|
        
        The following attributes are passed from ``self.attr`` to the new |SegmentChain|
        
            #. transcript_id, taken from :py:meth:`SegmentChain.get_name`
            #. gene_id, taken from :py:meth:`SegmentChain.get_gene`
            #. ID, generated as `"%s_3UTR" % self.get_name()`


        Parameters
        ----------
        extra_attr : keyword arguments
            Values that will be included in the 3' UTR subchain's `attr` dict.
            These can be used to overwrite values already present.


        Returns
        -------
        |SegmentChain|
            3' UTR region of `self` if present, otherwise empty |SegmentChain|
        """
        cdef:
            SegmentChain my_segmentchain
            dict attr = {}

        if self.cds_genome_start is not None and self.cds_genome_end is not None:

            my_segmentchain = self.c_get_subchain(self.cds_end,self.length,True)
            attr["type"] = "3UTR"
            attr["gene_id"] = self.get_gene()
            attr["transcript_id"] = self.get_gene()
            attr["ID"] = "%s_3UTR" % self.get_name()
            my_segmentchain.attr.update(attr)
            my_segmentchain.attr.update(extra_attr)

            return my_segmentchain  
        else:
            return SegmentChain()

    def as_gtf(self, str feature_type="exon", bint escape=True, list excludes=None):
        """Format `self` as a `GTF2`_ block. |GenomicSegments| are formatted
        as `GTF2`_ `'exon'` features. Coding regions, if peresent, are formatted
        as `GTF2`_ `'CDS'` features. Stop codons are excluded in the `'CDS'` features,
        per the `GTF2`_ specification, and exported separately.

        All attributes from `self.attr` are propagated to the exon and CDS
        features that are generated.

         
        Parameters
        ----------
        feature_type : str
            If not None, overrides the `'type'` attribute of `self.attr`
        
        escape : bool, optional
            URL escape tokens in column 9 of `GTF2`_ output (Default: `True`)
        
        
        Returns
        -------
        str
            Block of GTF2-formatted text


        Notes
        -----
        `gene_id` and `transcript_id` are required
            The `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_ requires
            that attributes `gene_id` and `transcript_id` be defined. If these
            are not present in `self.attr`, their values will be guessed 
            following the rules in :py:meth:`SegmentChain.get_gene` and 
            :py:meth:`SegmentChain.get_name`, respectively.
        
        Beware of attribute loss
            To save memory, only the attributes shared by all of the individual
            sub-features (e.g. exons) that were used to assemble this |Transcript|
            have been stored in `self.attr`. This means that upon re-export to `GTF2`_,
            these sub-features will be lacking any attributes that were specific
            to them individually. Formally, this is compliant with the 
            `GTF2 specification <http://mblab.wustl.edu/GTF22.html>`_, which states
            explicitly that only the attributes `gene_id` and `transcript_id`
            are supported.
            
        Columns of `GTF2`_ are as follows:
        
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes. "gene_id" and "transcript_id" are required                        
            ======== =========
        
        For more info
            - `GTF2 file format specification <http://mblab.wustl.edu/GTF22.html>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_        
        """
        cdef:
            str stmp
            SegmentChain cds_chain_temp, start_codon_chain, stop_codon_chain
            dict child_chain_attr
            list cds_positions
            GenomicSegment span = self.spanning_segment
            Strand my_strand = span.c_strand
            my_chrom = span.chrom

        excludes = [] if excludes is None else excludes

        stmp  = SegmentChain.as_gtf(self,feature_type=feature_type,escape=escape,excludes=[])
        cds_chain_temp = self.get_cds()
        if len(cds_chain_temp) > 0:
            child_chain_attr  = copy.deepcopy(self.attr)
            child_chain_attr.pop("type")
            cds_positions = list(cds_chain_temp.get_position_list())
            
            # remove stop codons from CDS, per GTF2 spec
            if my_strand == forward_strand:
                cds_positions = cds_positions[:-3]
            else:
                cds_positions = cds_positions[3:]
            cds_chain = SegmentChain(*positionlist_to_segments(my_chrom,
                                                               strand_to_str(my_strand),
                                                               cds_positions),
                                   type="CDS",**child_chain_attr)
            stmp += cds_chain.as_gtf(feature_type="CDS",escape=escape,excludes=excludes)
            
            start_codon_chain = cds_chain_temp.get_subchain(0, 3)
            start_codon_chain.attr.update(child_chain_attr)
            stmp += start_codon_chain.as_gtf(feature_type="start_codon",escape=escape,excludes=excludes)
    
            tmplen = cds_chain_temp.length
            stop_codon_chain = cds_chain_temp.get_subchain(tmplen-3, tmplen)
            stop_codon_chain.attr.update(child_chain_attr)
            stmp += stop_codon_chain.as_gtf(feature_type="stop_codon",escape=escape,excludes=excludes)

        return stmp

    def as_gff3(self, bint escape=True, list excludes=None, str rna_type="mRNA"):
        """Format a |Transcript| as a block of `GFF3`_ output, following
        the schema set out in the `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_
        
        The |Transcript| will be formatted according to the following rules:
        
          1. A feature of type `rna_type` will be created, with `Parent` attribute
             set to the value of ``self.get_gene()``, and `ID` attribute
             set to ``self.get_name()``
        
          2. For each |GenomicSegment| in `self`, a child feature of type
             `exon` will be created. The `Parent` attribute of these features
             will be set to the value of ``self.get_name()``. These will
             have unique IDs generated from ``self.get_name()``.

          3. If `self` is coding (i.e. has none-`None` value for
             `self.cds_genome_start` and `self.cds_genome_end`), child features
             of type `'five_prime_UTR'`, `'CDS'`, and `'three_prime_UTR'` will be created,
             with `Parent` attributes set to ``self.get_name()``. These will
             have unique IDs generated from ``self.get_name()``.
        
        
        Parameters
        ----------
        escape : bool, optional
            Escape tokens in column 9 of `GFF3`_ output (Default: `True`)
        
        excludes : list, optional
            List of attribute key names to exclude from column 9
            (Default: `[]`)
        
        rna_type : str, optional
            Feature type to export RNA as (e.g. `'tRNA'`, `'noncoding_RNA'`,
            et c. Default: `'mRNA'`)

        
        Returns
        -------
        str
            Multiline block of `GFF3`_-formatted text


        Notes
        -----
        Beware of attribute loss
            This |Transcript| was assembled from multiple individual component
            features (e.g. single exons), which may or may not have had their own
            unique attributes in their original annotation. To reduce overhead, 
            these individual attributes (if they were present) have not been
            (entirely) stored, and consequently will not (all) be exported.
            If this poses problems, consider instead importing, modifying, and
            exporting the component features

        GFF3 schemas vary
            Different GFF3s have different schemas (parent-child relationships
            between features). Here we adopt the commonly-used schema set by
            `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_,
            which may or may not match your schema.

        Columns of `GFF3`_ are as follows
            ======== =========
            Column   Contains
            ======== =========
                1     Contig or chromosome 
                2     Source of annotation 
                3     Type of feature ("exon", "CDS", "start_codon", "stop_codon") 
                4     Start (1-indexed)  
                5     End (fully-closed)
                6     Score  
                7     Strand  
                8     Frame. Number of bases within feature before first in-frame codon (if coding) 
                9     Attributes                       
            ======== =========

        For futher information, see
            - `GFF3 file format specification <http://www.sequenceontology.org/gff3.shtml>`_
            - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`
            - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_
            - `UCSC file format FAQ <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        cdef:
            str gene_id       = self.get_gene()
            str transcript_id = self.get_name()
            str ftype, my_id
            list ltmp = []
            list parts = []
            dict child_attr = copy.deepcopy(self.attr)
            SegmentChain feature
            GenomicSegment iv
            int n

        excludes = [] if excludes is None else excludes

        keys_to_pop = ("ID",)
        for k in keys_to_pop:
            if k in child_attr:
                child_attr.pop(k)

        # mRNA feature
        feature = SegmentChain(self.spanning_segment,ID=transcript_id,Parent=gene_id,type=rna_type)
        ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        # child features
        child_attr["Parent"] = transcript_id

        # exon feature
        child_attr["type"] = "exon"
        for n, seg in enumerate(self._segments):
            my_id   = "%s:exon:%s" % (transcript_id,n)
            child_attr["ID"]   = my_id
            feature = SegmentChain(seg,**child_attr)
            ltmp.append(SegmentChain.as_gff3(feature,excludes=excludes,escape=escape))
        
        # CDS & UTRs
        if self.cds_genome_start is not None:
            parts = [("five_prime_UTR", self.get_utr5()),
                     ("CDS",            self.get_cds()),
                     ("three_prime_UTR",self.get_utr3()),
                    ]
            for ftype, feature in parts:
                child_attr["type"]   = ftype
                child_attr["Parent"] = transcript_id
                for n, seg in enumerate(feature):
                    my_id   = "%s:%s:%s" % (transcript_id,ftype,n)
                    child_attr["ID"]   = my_id
                    feature = SegmentChain(seg,**child_attr)
                    ltmp.append(feature.as_gff3(excludes=excludes,escape=escape))

        return "".join(ltmp)
    
    def as_bed(self,as_int=True,color=None,extra_columns=None,empty_value=""):
        """Format `self` as a BED12[+X] line, assigning CDS boundaries 
        to the thickstart and thickend columns from `self.attr`

        If the |SegmentChain| was imported as a `BED`_ file with extra columns,
        these will be output in the same order, after the `BED`_ columns.
        
        Parameters
        ----------
        as_int : bool, optional
            Force "score" to integer (Default: True)
    
        color : str or None, optional
            Color represented as RGB hex string.
            If not none, overrides the color in `self.attr["color"]`
    
        extra_columns : None or list-like, optional
            If `None`, and the |SegmentChain| was imported using the `extra_columns`
            keyword of :meth:`~plastid.genomics.roitools.SegmentChain.from_bed`,
            the |SegmentChain| will be exported in BED 12+X format, in which
            extra columns are in the same order as they were upon import. If no extra columns
            were present, the |SegmentChain| will be exported as a BED12 line.

            If a list of attribute names, these attributes will be exported as
            extra columns in order, overriding whatever happened upon import. 
            If an attribute name is not in the `attr` dict of the |SegmentChain|,
            it will be exported with the value of `empty_value`

            If an empty list, no extra columns will be exported; the |SegmentChain|
            will be formatted as a BED12 line.
        
        empty_value : str, optional

    
        Returns
        -------
        str
            Line of BED12-formatted text
            
        
        Notes
        -----
        BED12 columns are as follows
            ======== =========
            Column   Contains
            ======== =========
               0     Contig or chromosome
               1     Start of first block in feature (0-indexed)
               2     End of last block in feature (half-open)
               3     Feature name
               4     Feature score
               5     Strand
               6     thickstart
               7     thickend
               8     Feature color as RGB tuple
               9     Number of blocks in feature
               10    Block lengths
               11    Block starts, relative to start of first block
            ======== =========

        Fore more information
            See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
        """
        return SegmentChain.as_bed(self,
                                   thickstart=self.cds_genome_start,
                                   thickend=self.cds_genome_end,
                                   as_int=as_int,
                                   color=color,
                                   extra_columns=extra_columns,
                                   empty_value=empty_value)

    @staticmethod
    def from_bed(str line, extra_columns=0):
        """Create a |Transcript| from a BED line with 4 or more columns.
        `thickstart` and `thickend` columns, if present, are assumed to specify
        CDS boundaries, a convention that, while common, is formally outside the
        `BED`_ specification.
    
    	See the `UCSC file format faq <http://genome.ucsc.edu/FAQ/FAQformat.html>`_
    	for more details.

        Parameters
        ----------
        line
            Line from a BED file with at least 4 columns

        extra_columns: int or list, optional
            Extra, non-BED columns in `BED`_ file corresponding to feature
            attributes. This is common in `ENCODE`_-specific `BED`_ variants.
            
            if `extra-columns` is:
            
              - an :class:`int`: it is taken to be the
                number of attribute columns. Attributes will be stored in
                the `attr` dictionary of the |SegmentChain|, under names like
                `custom0`, `custom1`, ... , `customN`.

              - a :class:`list` of :class:`str`, it is taken to be the names
                of the attribute columns, in order, from left to right in the file.
                In this case, attributes in extra columns will be stored under

              - a :class:`list` of :class:`tuple`, each tuple is taken
                to be a pair of `(attribute_name, formatter_func)`. In this case,
                the value of `attribute_name` in the `attr` dict of the |SegmentChain|
                will be set to `formatter_func(column_value)`.
            
            (Default: 0)
                
    
        Returns
        -------
        |Transcript|
        """
        cdef:
            dict         attr     = get_attr_from_bed(line,extra_columns=extra_columns)
            list         segments = get_segments_from_bed(attr)
            Transcript   chain    = Transcript()

        chain._set_segments(segments)
        cds_genome_start = attr.get("thickstart",None)
        cds_genome_end   = attr.get("thickend",None)
 
        if cds_genome_start == cds_genome_end:
            cds_genome_start = cds_genome_end = None
 
        chain.cds_genome_start = cds_genome_start
        chain.cds_genome_end   = cds_genome_end
 
        if chain.cds_genome_start is not None and chain.cds_genome_end is not None:
            chain._update_cds()
 
        attr["type"] = "mRNA" # default type for SegmentChain is "exon". We want to use "mRNA"
        attr.pop("thickstart")
        attr.pop("thickend")

        chain.attr = attr
        
        return chain
    
    @staticmethod
    def from_psl(str psl_line):
        cdef:
            Transcript transcript = Transcript()
            SegmentChain segchain = SegmentChain.from_psl(psl_line)
            dict attr = segchain.attr

        attr["cds_genome_start"] = None
        attr["cds_genome_end"]   = None
        transcript._set_segments(segchain._segments)
        transcript.attr = attr

        return transcript


