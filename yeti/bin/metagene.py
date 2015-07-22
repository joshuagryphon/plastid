#!/usr/bin/env python
"""Performs :term:`metagene` analyses. The workflow is separated into the
following subprograms:

Generate
    A :term:`metagene` profile is a position-wise average over all genes
    in the vicinity of an interesting landmark (e.g. a start codon). Because
    genes can have multiple transcript isoforms that may cover different
    genomic positions, which transcript positions (and therefore which
    genomic positions) to include in the average can be ambiguous.
    
    To avoid this problem, we define for each gene a *maximum spanning window*
    of unambiguous positions surrounding the landmark of interest, and save
    these windows into an ROI file. The windows are defined by the following
    algorithm: 
    
    1.  Transcripts are grouped by gene. If all transcripts share the same
        genomic coordinate for the landmark of interest (e.g. if all 
        transcripts share the same start codon), then all transcripts are
        included in the analysis. If not, all transcripts and their associated
        gene are excluded from further processing.
    
    2.  For each set of transcripts that pass step (1), the maximal spanning
        window is created by aligning the set of transcripts at the landmark, and
        adding nucleotide positions in transcript coordinates to the growing
        window in both 5' and 3' directions until either:
        
            - the next nucleotide position added is no longer corresponds to 
              the same genomic position in all transcripts
            
            - the window reaches the maximum user-specified size

    **Note**: if annotations are supplied as `BED`_ files, transcripts cannot be
    grouped by gene, because `BED`_ files don't contain this information.
    In this case one ROI is generated per transcript. This may or may not
    be what you want. You can filter later and decide.
    
    
    .. Rubric :: Output files
    
    OUTBASE_rois.txt
        A tab-delimited text file describing the maximal spanning window for
        each gene, with columns as follows:
        
        ================   ==================================================
        **Column**         **Contains**
        ----------------   --------------------------------------------------

        alignment_offset   Offset to align window to all other windows in the
                           file, if the window happens to be shorter on the 5'
                           end than specified in ``--flank_upstream``. Typically
                           this is `0`.

        gene_id            ID of gene
        
        region             maximal spanning window, formatted as
                           `chromosome:start-end:(strand)`
        
        window_size        with of window
        
        zero_point         distance from 5' end of window to landmark
        ================   ==================================================
        
    
    OUTBASE_rois.bed
        Maximal spanning windows in `BED`_ format for visualization in
        a :term:`genome browser`. The thickly-rendered portion of a window
        indicates its landmark

    where `OUTBASE` is supplied by the user.
    
    
Count
    This program generate :term:`metagene` profiles from :term:`counts` or
    :term:`alignments`, taking the following steps:
    
    1.  The **raw counts** at each position in each window (from the ``generate``
        subprogram) are totaled to create a raw count vector for the window.

    2.  A **normalized count vector** is created for each window by dividing
        its raw count vector by the total number of counts occurring within a
        user-defined normalization window within the window.
    
    3.  A **metagene average** is created by taking aligning all of the
        normalized count vectors, and taking the median normalized counts
        over all vectors at each nucleotide position. Count vectors deriving
        from windows that don't meet a minimum count threshold (set via the
        ``--norm_region`` option) are excluded.
    
    
    .. Rubric :: Output files

    Raw count vectors, normalized count vectors, and metagene profiles are all
    saved as tab-delimited text files, for subsequent plotting, filtering,
    or reanalysis.
    
    OUTBASE_metagene_profile.txt
        Tab-delimited table of metagene profile, containing the following
        columns:

        ================   ==================================================
        **Column**         **Contains**
        ----------------   --------------------------------------------------
        x                  Distance in nucleotides from the landmark
        
        metagene_average   Value of metagene average at that position
        
        regions_counted    Number of maximal spanning windows included at
                           that position in the average. i.e. windows that
                           both met the threshold set by ``--min_counts`` and
                           were not masked at that position by a :term:`mask file`
        ================   ==================================================        
        
    OUTBASE_rawcounts.txt
        Table of raw counts. Each row is a maximal spanning window for a gene,
        and each column a nucleotide position in that window. All windows
        are aligned at the landmark.
    
    OUTBASE_normcounts.txt
        Table of normalized counts, produced by dividing each row in the
        raw count table by the of counts in that row within the columns
        specified by ``--norm_region``.
    
    where `OUTBASE` is supplied by the user.

    
Chart
    One or more metagene profiles generated by the ``count`` subprogram,
    for example, on different datasets, are plotted against each other. 


See command-line help for each subprogram for details on parameters for each 
"""


__author__ = "joshua"

import sys
import argparse
import numpy
from yeti.genomics.roitools import SegmentChain, positionlist_to_segments
from yeti.util.array_table import ArrayTable
from yeti.util.io.filters import NameDateWriter
from yeti.util.io.openers import get_short_name, argsopener, NullWriter
from yeti.util.scriptlib.argparsers import get_genome_array_from_args,\
                                                      get_transcripts_from_args,\
                                                      get_alignment_file_parser,\
                                                      get_annotation_file_parser,\
                                                      get_mask_file_parser,\
                                                      get_genome_hash_from_mask_args
from yeti.util.scriptlib.help_formatters import format_module_docstring
from yeti.util.services.decorators import catch_warnings

import inspect
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))

#===============================================================================
# helper functions to generate/handle maximum spanning windows / ROIs
#===============================================================================

def window_landmark(region,flank_upstream=50,flank_downstream=50,ref_delta=0,landmark=0):
    """Define a window surrounding a landmark in a region, if the region has such a landmark,
    (e.g. a start codon in a transcript), accounting for splicing of the region,
    if the region is discontinuous
    
    Parameters
    ----------
    transcript : |SegmentChain| or |Transcript|
        Region on which to generate a window surrounding a landmark
    
    landmark : int
        Position of the landmark within `region`

    flank_upstream : int
        Nucleotides upstream of `landmark` to include in window
    
    flank_downstream : int
        Nucleotides downstream of `landmark` to include in window
    
    ref_delta : int
        Offset from `landmark` to the reference point. If 0, the landmark
        is the reference point. Default: 0
    
    
    Returns
    -------
    |SegmentChain|
        Window of `region` surrounding landmark
    
    int
        alignment offset to the window start, if `region` itself wasn't long
        enough in the 5' direction to include the entire distance specified by
        `flank_upstream`. Use this to align windows generated around similar
        landmarks from different `regions` (e.g. windows surrounding start codons
        in various transcripts).

    (str, int, str)
        Genomic coordinate of reference point as *(chromosome name, coordinate, strand)*
    """
    if landmark + ref_delta >= flank_upstream:
        fiveprime_offset   = 0
        my_start = landmark + ref_delta - flank_upstream
    else:
        fiveprime_offset  = flank_upstream - landmark
        my_start = 0
    
    my_end = min(region.get_length(),landmark + ref_delta + flank_downstream)
    roi    = region.get_subchain(my_start,my_end)
   
    if landmark + ref_delta == region.get_length():
        if region.spanning_segment.strand == "+":
            ref_point = (region.spanning_segment.chrom,region.spanning_segment.end,region.spanning_segment.strand)
        else:
            ref_point = (region.spanning_segment.chrom,region.spanning_segment.start - 1,region.spanning_segment.strand)
    else:
        ref_point = region.get_genomic_coordinate(landmark+ref_delta)

    return roi, fiveprime_offset, ref_point
    
def window_cds_start(transcript,flank_upstream,flank_downstream,ref_delta=0):
    """Returns a window surrounding a start codon.
    
    Parameters
    ----------
    transcript : |Transcript|
        Transcript on which to generate window
    
    flank_upstream : int
        Nucleotide length upstream of start codon to include in window,
        if `transcript` has a start codon
    
    flank_downstream : int
        Nucleotide length downstream of start codon to include in window,
        if `transcript` has a start codon
    
    ref_delta : int, optional
        Offset from  start codon to the reference point. If `0`, the landmark
        is the reference point. (Default: `0`)
    
    Returns
    -------
    |SegmentChain|
        Window surrounding start codon if `transcript` is coding. Otherwise,
        zero-length |SegmentChain| 
    
    int
        alignment offset to the window start, if `transcript` itself wasn't long
        enough in the 5\' direction to include the entire distance specified by
        `flank_upstream`. Use this to align this window to other windows generated
        around start codons in other transcripts.

    (str, int, str)
        Genomic coordinate of reference point as *(chromosome name, coordinate, strand)*
    """
    if transcript.cds_start is None:
        return SegmentChain(), numpy.nan, numpy.nan

    return window_landmark(transcript,flank_upstream,flank_downstream,
                           ref_delta=ref_delta,
                           landmark=transcript.cds_start)

def window_cds_stop(transcript,flank_upstream,flank_downstream,ref_delta=0):
    """Returns a window surrounding a stop codon.

    Parameters
    ----------
    transcript : |Transcript|
        Transcript on which to generate window
    
    flank_upstream : int
        Nucleotide length upstream of stop codon to include in window,
        if `transcript` has a stop codon
    
    flank_downstream : int
        Nucleotide length downstream of stop codon to include in window,
        if `transcript` has a stop codon
    
    ref_delta : int, optional
        Offset from  stop codon to the reference point. If `0`, the landmark
        is the reference point. (Default: `0`)
    
    Returns
    -------
    |SegmentChain|
        Window surrounding stop codon if transcript is coding. Otherwise,
        zero-length |SegmentChain| 
    
    int
        alignment offset to the window start, if `transcript` itself wasn't long
        enough in the 5' direction to include the entire distance specified by
        `flank_upstream`. Use this to align this window to other windows generated
        around stop codons in other transcripts.

    (str, int, str)
        Genomic coordinate of reference point as *(chromosome name, coordinate, strand)*
    """
    if transcript.cds_start is None:
        return SegmentChain(), numpy.nan, numpy.nan
    
    return window_landmark(transcript,flank_upstream,flank_downstream,
                           ref_delta=ref_delta,
                           landmark=transcript.cds_end-3)


#===============================================================================
# Subprograms
#===============================================================================

def do_generate(transcripts,mask_hash,flank_upstream,flank_downstream,
                landmark_func=window_cds_start,
                printer=NullWriter()):
    """Generate a file of maximal spanning windows surrounding a landmark,
for use in ``count`` subprogram. Windows are generated by the following 
algorithm:

    1.  Transcripts are grouped by gene. If all transcripts share the same
        genomic coordinate for the landmark of interest (e.g. if all 
        transcripts share the same start codon), then all transcripts are
        included in the analysis. If not, all transcripts and their associated
        gene are excluded from further processing.
    
    2.  For each set of transcripts that pass step (1), the maximal spanning
        window is created by aligning the set of transcripts at the landmark, and
        adding nucleotide positions in transcript coordinates to the growing
        window in both 5\' and 3\' directions until either:
        
            - the next nucleotide position added is no longer corresponds to 
              the same genomic position in all transcripts
            
            - the window reaches the maximum user-specified size
    
    Parameters
    ----------
    transcripts : list
        A list of |Transcript| objects, preferably with `gene_id` and
        `transcript_id` (e.g. transcripts assembled from a `GTF2`_ or `GFF3`_
        file), so that transcripts can be grouped by gene when making maximal
        spanning windows.
    
    mask_hash : |GenomeHash|
        |GenomeHash| containing regions to exclude from analysis
    
    flank_upstream : int
        Number of nucleotides upstream of landmark to include in windows
        (in transcript coordinates)

    flank_downstream: int
        Number of nucleotides downstream of landmark to include in windows
        (in transcript coordinates)
    
    landmark_func : func, optional
        Function yielding coordinate of landmark in transcript. As examples,
        :py:func:`window_cds_start` and :py:func:`window_cds_stop` are provided,
        though any function that meets the following criteria can be used:
        
            1. It must take the same parameters as :func:`window_cds_start`
            
            2. It must return the same types as :func:`window_cds_start`

        Such functions could choose arbitrary features as landmarks, such as
        peaks in ribosome density, nucleic acid sequence features, transcript
        start or end sites, or any property that can be deduced from a
        |Transcript|. (Default: :func:`window_cds_start`)
    
    printer : file-like, optional
        filehandle to write logging info to (Default: :func:`NullWriter`)
    
    Returns
    -------
    |ArrayTable|
        An |ArrayTable| containing the following columns describing the
        maximal spanning windows:

            ================   ==================================================
            *Column*           *Contains*
            ----------------   --------------------------------------------------
    
            alignment_offset   Offset to align window to all other windows in the
                               file, if the window happens to be shorter on the 5\'
                               end than specified in ``--flank_upstream``
    
            gene_id            ID of gene
            
            region             maximal spanning window, formatted as
                               `chromosome:start-end:(strand)`
            
            window_size        with of window
            
            zero_point         distance from 5' end of window to landmark
            ================   ==================================================        
    
    list
        List of |SegmentChain| representing each window. These data are also
        represented as strings in the |ArrayTable|
    
    Notes
    -----         
    Not all genes will be included in the output if, for example, there isn't a
    position set common to all transcripts surrounding the landmark
    """
    window_size = flank_upstream + flank_downstream
        
    dtmp = { "gene_id"           : [],
             "region"            : [],
             "masked"            : [],
             "alignment_offset"  : [],
             "window_size"       : [],
             "zero_point"        : [],
             }
    
    # build gene-transcript graph
    printer.write("Building gene-transcript graph...")
    gene_transcript = {}
    for tx_chain in transcripts:
        gene_id = tx_chain.get_gene()
        try:
            gene_transcript[gene_id].append(tx_chain)
        except KeyError:
            gene_transcript[gene_id] = [tx_chain]
    
    # for each gene, find maximal window in which all points
    # are represented in all transcripts. return IVC and offset
    c = -1
    export_rois = []
    for gene_id, txlist in sorted(gene_transcript.items()):
        c += 1
        if c % 1000 == 1:
            printer.write("Processed %s transcripts, included %s..." % (c,len(list(dtmp.values())[0])))
        my_roi    = None
        new_roi   = None
        my_offset = None
        new_roi = zero_point_roi = None
        new_offset = None
        refpoints = []
        
        # find common positions
        position_matrix = numpy.tile(numpy.nan,(len(txlist),window_size))
        for n,tx in enumerate(txlist):
            try:
                my_roi, my_offset, genomic_refpoint = landmark_func(tx,flank_upstream,flank_downstream)
                refpoints.append(genomic_refpoint)
                
                if genomic_refpoint is not numpy.nan and len(my_roi) > 0:
                    pos_list = my_roi.get_position_list() # ascending list of positions
                    my_len = len(pos_list)
                    assert my_offset + my_len <= window_size 
                    if my_roi.spanning_segment.strand == "+":
                        position_matrix[n,my_offset:my_offset+my_len] = pos_list
                    else:
                        my_len = len(pos_list)
                        position_matrix[n,my_offset:my_offset+my_len] = pos_list[::-1]
                
            except IndexError:
                printer.write("IndexError at %s, transcript %s: " % (gene_id,tx.get_name()))
        
        # continue only if refpoints all match
        if len(set(refpoints)) == 1 and numpy.nan not in refpoints:
            new_shared_positions = []
            if len(set(refpoints)) == 1:
                for i in range(0,position_matrix.shape[1]):
                    col = position_matrix[:,i]
                    if len(set(col)) == 1 and col[0] is not numpy.nan:
                        new_shared_positions.append(col[0])
          
            # continue only if there exist positions shared between all regions 
            new_shared_positions = [int(X) for X in new_shared_positions if not numpy.isnan(X)]
            if len(set(new_shared_positions)) > 0:
        
                # define new ROI covering all positions common to all transcripts
                new_roi = SegmentChain(*positionlist_to_segments(txlist[0].chrom,
                                                                 txlist[0].strand,
                                                                 new_shared_positions))
    
                if new_roi.spanning_segment.strand == "+":
                    new_roi.attr["thickstart"] = genomic_refpoint[1]
                    new_roi.attr["thickend"]   = genomic_refpoint[1] + 1
                else:
                    new_roi.attr["thickstart"] = genomic_refpoint[1]
                    new_roi.attr["thickend"]   = genomic_refpoint[1] + 1
    
                # having made sure that refpoint is same for all transcripts,
                # we use last ROI and last offset to find new offset
                # this fails if ref point is at the 3' end of the roi, 
                # due to quirks of half-open coordinate systems
                # so we test it explicitly
                if flank_upstream - my_offset == my_roi.get_length():
                    new_offset = my_offset
                else:
                    zero_point_roi = new_roi.get_segmentchain_coordinate(*genomic_refpoint)
                    new_offset = flank_upstream - zero_point_roi
        
                masks = mask_hash.get_overlapping_features(new_roi)
                mask_segs = []
                for mask in masks:
                    mask_segs.extend(mask._segments)
                
                mask_chain = SegmentChain(*mask_segs)
                
                dtmp["gene_id"].append(gene_id)
                dtmp["window_size"].append(window_size)
                dtmp["region"].append(str(new_roi)) # need to cast to string to keep numpy from converting to array
                dtmp["masked"].append(str(mask_chain))
                dtmp["alignment_offset"].append(new_offset)
                dtmp["zero_point"].append(flank_upstream)
                export_rois.append(new_roi)

    # convert to ArrayTable
    dtmp = ArrayTable(dtmp)
    printer.write("Processed %s transcripts total. Included %s." % (c+1,len(dtmp)))

    return dtmp, export_rois

@catch_warnings("ignore")
def do_count(roi_table,ga,norm_start,norm_end,min_counts,printer=NullWriter()):
    """Calculate a metagene average over maximal spanning windows specified in 
    `roi_table`, taking the following steps:
    
    1.  The **raw counts** at each position in each window (from the ``generate``
        subprogram) are totaled to create a raw count vector for the ROI.

    2.  A **normalized count vector** is created fore each window by dividing
        its raw count vector by the total number of counts occurring within a
        user-defined normalization window within the window.
    
    3.  A **metagene average** is created by taking aligning all of the
        normalized count vectors, and taking the median normalized counts
        over all vectors at each nucleotide position. Count vectors deriving
        from ROIs that don't meet a minimum count threshold (set via the
        ``--norm_region`` option) are excluded.
            
    Parameters
    ----------
    roi_table : |ArrayTable|
        |ArrayTable| specifying maximal spanning windows, generated
        by :py:meth:`do_generate`
    
    ga : instance of subclass of |AbstractGenomeArray|
        Count or alignment data
    
    norm_start : int
        Coordinate in window specifying normalization region start
    
    norm_end : int
        Coordinate in window specifying normalization region end
    
    min_counts : float
        Minimum number of counts in `window[norm_start:norm_end]`
        required for inclusion in metagene profile

    printer : file-like, optional
        filehandle to write logging info to (Default: :func:`NullWriter`)

        
    Returns
    -------
    :py:class:`numpy.ndarray`
        raw counts at each position (column) in each window (row)
    
    :py:class:`numpy.ndarray`
        counts at each position (column) in each window (row), normalized by
        the total number of counts in that row from `norm_start` to `norm_end`
    
    |ArrayTable|
        Metagene profile of median normalized counts at each position across
        all windows, and the number of windows included in the calculation of each
        median
    """
    window_size    = roi_table["window_size"][0]
    upstream_flank = roi_table["zero_point"][0]
    counts = numpy.ma.MaskedArray(numpy.tile(numpy.nan,(len(roi_table),window_size)))
    
    for i in range(len(roi_table)):
        if i % 1000 == 1:
            printer.write("Counted %s ROIs..." % (i))
            
        roi    = SegmentChain.from_str(roi_table["region"][i])
        mask   = SegmentChain.from_str(roi_table["masked"][i])
        roi.add_masks(*mask)
        offset = int(round((roi_table["alignment_offset"][i])))
        assert offset + roi.get_length() <= window_size
        counts[i,offset:offset+roi.get_length()] = roi.get_valid_counts(ga)
    
    printer.write("Counted %s ROIs total." % (i+1))
        
    denominator = numpy.nansum(counts[:,norm_start:norm_end],1)
    norm_counts = (1.0*counts.T / denominator).T
    
    norm_counts = numpy.ma.masked_invalid(norm_counts)
    
    profile   = numpy.ma.median(norm_counts[denominator >= min_counts],axis=0)
    num_genes = ((~norm_counts.mask)[denominator >= min_counts]).sum(0) 
    
    profile_table = ArrayTable({ "metagene_average"     : profile,
                                 "regions_counted"      : num_genes,
                                 "x"                    : numpy.arange(-upstream_flank,window_size-upstream_flank),
                                })
    return counts, norm_counts, profile_table


def do_chart(sample_dict,landmark="landmark",title=None):
    """Plot metagene profiles
    
    Parameters
    ----------
    sample_dict : dict
        Dictionary mapping sample names to |ArrayTable| s containing metagene
        profile information from ``count`` subprogram
    
    landmark : str, optional
        Name of landmark at zero point, used for labeling X-axis. e.g.
        `'Start codon'` (Default: `'landmark'`)
        
    title : str or None, optional
        Chart title (Default: `None`)
    
    Returns
    -------
    :py:class:`matplotlib.Figure`
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure()
    min_x = numpy.inf
    max_x = -numpy.inf
    for k, v in sample_dict.items():
        plt.plot(v["x"],v["metagene_average"],label=k)
        min_x = min(min_x,min(v["x"]))
        max_x = max(max_x,max(v["x"]))
    
    plt.xlim(min_x,max_x)
    ylim = plt.gca().get_ylim()
    plt.ylim(0,ylim[1])
    
    plt.xlabel("Distance from %s (nt)" % landmark)
    plt.ylabel("Normalized read density (au)")
    
    if title is not None:
        plt.title(title)
    
    return fig
        

#===============================================================================
# PROGRAM BODY
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
    
    alignment_file_parser  = get_alignment_file_parser(disabled=["normalize"])
    annotation_file_parser = get_annotation_file_parser()
    mask_file_parser = get_mask_file_parser()
    
    generator_help = "Create ROI file from genome annotation"
    generator_desc = format_module_docstring(do_generate.__doc__)
    
    count_help = "Count reads falling into regions of interest, normalize, and average into a metagene profile"
    count_desc = format_module_docstring(do_count.__doc__)

    chart_help = "Plot metagene profiles"
    chart_desc = format_module_docstring(do_chart.__doc__)

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title="subcommands",
                                       description="choose one of the following",
                                       dest="program")
    gparser    = subparsers.add_parser("generate",
                                       help=generator_help,
                                       description=generator_desc,
                                       parents=[annotation_file_parser,mask_file_parser],
                                       formatter_class=argparse.RawDescriptionHelpFormatter)
    cparser    = subparsers.add_parser("count",
                                       help=count_help,
                                       description=count_desc,
                                       parents=[alignment_file_parser],
                                       formatter_class=argparse.RawDescriptionHelpFormatter)
    pparser    = subparsers.add_parser("chart",
                                       help=chart_help,
                                       description=chart_desc,
                                       formatter_class=argparse.RawDescriptionHelpFormatter)
    
    # generate subprogram options
    gparser.add_argument("--landmark",type=str,choices=("cds_start","cds_stop"),
                         default="cds_start",
                         help="Landmark around which to build metagene profile (Default: cds_start)")
    gparser.add_argument("--upstream",type=int,default=50,
                         help="Nucleotides to include upstream of landmark (Default: 50)")
    gparser.add_argument("--downstream",type=int,default=50,
                         help="Nucleotides to include downstream of landmark (Default: 50)")
    gparser.add_argument("outbase",type=str,
                         help="Basename for output files")
    
    # count subprogram options
    cparser.add_argument("roi_file",type=str,
                         help="Text file containing maximal spanning windows and offsets, "+
                              "generated by the ``generate`` subprogram.")
    cparser.add_argument("--min_counts",type=int,default=10,metavar="N",
                         help="Minimum counts required in normalization region "+
                              "to be included in metagene average (Default: 10)")
    cparser.add_argument("--norm_region",type=int,nargs=2,metavar="N",
                         default=(70,100),
                         help="Portion of each window against which its individual raw count profile"+
                              " will be normalized. Specify two integers, in nucleotide"+
                              " distance, from 5\' end of window. (Default: 70 100)")
    cparser.add_argument("outbase",type=str,
                         help="Basename for output files")

    # chart subprogram arguments
    pparser.add_argument("outfile",type=str,
                         help="Name of output file. Format will be auto-detected"+
                         " from the extension. Valid options are `png`, `svg`, `pdf`, `eps`, et c.")
    pparser.add_argument("infiles",type=str,nargs="+",
                         help="One or more metagene profiles, generated by the"+
                              " ``count`` subprogram, which will be plotted together."
                         )
    pparser.add_argument("--labels",type=str,nargs="+",default=[],
                         help="Sample names for each metagene profile (optional).")
    pparser.add_argument("--title",type=str,default=None,
                         help="Title for chart (optional)")
    pparser.add_argument("--landmark",type=str,default=None,
                         help="Name of landmark at zero point (e.g. 'CDS start' or 'CDS stop'; optional)")

    
    args = parser.parse_args(argv)
    
    # 'generate' subprogram
    if args.program == "generate":
        printer.write("Generating ROI file...")
        if args.landmark == "cds_start":
            map_function = window_cds_start
        elif args.landmark == "cds_stop":
            map_function = window_cds_stop
            
        # open annotations
        printer.write("Opening annotation files: %s..." % ", ".join(args.annotation_files))
        transcripts = get_transcripts_from_args(args,printer=printer)
        
        mask_hash = get_genome_hash_from_mask_args(args)
        
        # get ROIs
        printer.write("Generating regions of interest...")   
        roi_table, export_rois = do_generate(transcripts,
                                             mask_hash,
                                             args.upstream,
                                             args.downstream,
                                             landmark_func=map_function,
                                             printer=printer)
        roi_file = "%s_rois.txt" % args.outbase
        bed_file = "%s_rois.bed" % args.outbase
        printer.write("Saving to ROIs %s..." % roi_file)
        with argsopener(roi_file,args,"w") as roi_fh:
            roi_table.to_file(roi_fh)
            roi_fh.close()
            
        printer.write("Saving BED output as %s..." % bed_file)
        with argsopener(bed_file,args,"w") as bed_fh:
            for roi in export_rois:
                bed_fh.write(roi.as_bed())
        
            bed_fh.close()
    
    # 'count' subprogram
    elif args.program == "count":
        printer.write("Opening ROI file %s..." % args.roi_file)
        with open(args.roi_file) as fh:
            roi_table = ArrayTable.from_file(fh)
            fh.close()
        
        # open count files
        gnd = get_genome_array_from_args(args,printer=printer,disabled=["normalize"])
        
        # count
        printer.write("Counting...")
        counts, norm_counts, profile_table = do_count(roi_table,
                                                      gnd,
                                                      args.norm_region[0],
                                                      args.norm_region[1],
                                                      args.min_counts,
                                                      printer=printer)
        
        # save
        count_fn     = "%s_rawcounts.txt" % args.outbase
        normcount_fn = "%s_normcounts.txt" % args.outbase
        profile_fn   = "%s_metagene_profile.txt" % args.outbase
        printer.write("Saving counts to %s..." % count_fn)
        numpy.savetxt(count_fn,counts,delimiter="\t",fmt='%.8f')
        printer.write("Saving normalized counts to %s..." % normcount_fn)
        numpy.savetxt(normcount_fn,norm_counts,delimiter="\t")
        printer.write("Saving metagene profile to %s..." % profile_fn)
        with argsopener(profile_fn,args,"w") as profile_out:
            profile_table.to_file(profile_out,
                                  keyorder=["x","metagene_average","regions_counted"],
                                  formatters={numpy.float64  : '{:.8f}'.format,
                                              numpy.float32  : '{:.8f}'.format,
                                              numpy.float128 : '{:.8f}'.format,
                                  }
                                  )
            profile_out.close()
        
    # 'plot' subprogram
    elif args.program == "chart":
        assert len(args.labels) in (0,len(args.infiles))
        if len(args.labels) == len(args.infiles):
            samples = { K : ArrayTable.from_file(open(V)) for K,V in zip(args.labels,args.infiles)}
        else:
            samples = { get_short_name(V) : ArrayTable.from_file(open(V)) for V in args.infiles }
        
        figure = do_chart(samples,title=args.title,landmark=args.landmark)
        printer.write("Saving to %s..." % args.outfile)
        figure.savefig(args.outfile)

    printer.write("Done.")

if __name__ == "__main__":
    main()
