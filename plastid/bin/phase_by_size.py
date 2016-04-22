#!/usr/bin/env python
"""Estimate :term:`sub-codon phasing` in a :term:`ribosome profiling` dataset,
stratified by read length.

Because ribosomes step three nucleotides in each cycle of translation elongation,
in many :term:`ribosome profiling` datasets a triplet periodicity is observable
in the distribution of :term:`ribosome-protected footprints <footprint>`.

In a good dataset, 70-90% of the reads on a codon fall within the first of the
three codon positions. This allows deduction of translation reading frames, if
the reading frame is not known *a priori.* See :cite:`Ingolia2009` for more
details.

Output files
------------
    OUTBASE_phasing.txt
        Read phasing for each read length

    OUTBASE_phasing.svg
        Plot of phasing by read length

where `OUTBASE` is supplied by the user.

 .. note::

    To avoid double-counting of codons, users should either use an *ROI file*
    made by the ``generate`` subprogram of the |metagene| script, or supply
    an :term:`annotation` file that includes only one transcript isoform per
    gene.

"""
import sys
import argparse
import inspect
import warnings

import pandas as pd
import numpy
import matplotlib
matplotlib.use("Agg")
from plastid.util.scriptlib.argparsers import (AlignmentParser, AnnotationParser,
                                               PlottingParser, BaseParser)
from plastid.util.io.openers import get_short_name, argsopener, read_pl_table
from plastid.util.io.filters import NameDateWriter
from plastid.util.scriptlib.help_formatters import format_module_docstring
from plastid.util.services.exceptions import DataWarning, ArgumentWarning
from plastid.plotting.plots import phase_plot
from plastid.genomics.roitools import SegmentChain

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


def roi_row_to_cds(row):
    """Helper function to extract coding portions from maximal spanning windows
    flanking CDS starts that are created by |metagene| ``generate`` subprogram.
    
    Parameters
    ----------
    row : (int, Series)
        Row from a :class:`pandas.DataFrame` of an ROI file made by the |metagene|
        ``generate`` subprogram
        
    Returns
    -------
    |SegmentChain|
        Coding portion of maximal spanning window
    """
    chainstr, alignment_offset, zero_point = row[1][["region","alignment_offset","zero_point"]]
    chain = SegmentChain.from_str(chainstr)
    cds_start = zero_point - alignment_offset
    subchain = chain.get_subchain(cds_start,chain.length)
    return subchain

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
    al = AlignmentParser(disabled=["normalize","big_genome","spliced_bowtie_files"],
                        input_choices=["BAM"])
    an = AnnotationParser()
    pp = PlottingParser()
    bp = BaseParser()
    plotting_parser = pp.get_parser()
    
    alignment_file_parser = al.get_parser(conflict_handler="resolve")
    annotation_file_parser = an.get_parser(conflict_handler="resolve")
    base_parser = bp.get_parser()

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     conflict_handler="resolve",
                                     parents=[base_parser,
                                              annotation_file_parser,
                                              alignment_file_parser,
                                              plotting_parser])
    
    parser.add_argument("roi_file",type=str,nargs="?",default=None,
                        help="Optional. ROI file of maximal spanning windows surrounding start codons, "+\
                             "from ``metagene generate`` subprogram. Using this instead of `--annotation_files` "+\
                             "prevents double-counting of codons when multiple transcript isoforms exist "+\
                             "for a gene. See the documentation for `metagene` for more info about ROI files."+\
                             "If an ROI file is not given, supply an annotation with ``--annotation_files``")
    parser.add_argument("outbase",type=str,help="Required. Basename for output files")
    parser.add_argument("--codon_buffer",type=int,default=5,
                        help="Codons before and after start codon to ignore (Default: 5)")


    args = parser.parse_args(argv)
    bp.get_base_ops_from_args(args)
    pp.set_style_from_args(args)
    gnd = al.get_genome_array_from_args(args,printer=printer)

    read_lengths = list(range(args.min_length,args.max_length+1))
    codon_buffer = args.codon_buffer
    
    dtmp = { "read_length"   : numpy.array(read_lengths),
             "reads_counted" : numpy.zeros_like(read_lengths,dtype=int),
            }

    if args.roi_file is not None:
        using_roi = True
        roi_table = read_pl_table(args.roi_file)
        regions = roi_table.iterrows()
        transform_fn = roi_row_to_cds
        back_buffer = -1
        if len(args.annotation_files) > 0:
            warnings.warn("If an ROI file is given, annotation files are ignored. Pulling regions from '%s'. Ignoring '%s'" % (args.roi_file,
                                                                                                                               ", ".join(args.annotation_files)),
                          ArgumentWarning)
    else:
        using_roi = False
        if len(args.annotation_files) == 0:
            printer.write("Either an ROI file or at least annotation file must be given.")
            sys.exit(1)
        else:
            warnings.warn("Using a transcript annotation file instead of an ROI file can lead to double-counting of codons if the annotation contains multiple transcripts per gene.",
                          ArgumentWarning)        
            regions = an.get_transcripts_from_args(args,printer=printer)
            back_buffer  = -codon_buffer
            transform_fn = lambda x: x.get_cds()
    
    phase_sums = {}
    for k in read_lengths:
        phase_sums[k] = numpy.zeros(3)
    
    for n, roi in enumerate(regions):
        if n % 1000 == 1:
            printer.write("Counted %s ROIs ..." % n)
        
        # transformation needed to extract CDS from transcript or from ROI file window
        cds_part = transform_fn(roi)
        
        # only calculate for coding genes
        if len(cds_part) > 0:

            read_dict     = {}
            count_vectors = {}
            for k in read_lengths:
                read_dict[k]     = []
                count_vectors[k] = []
            
            # for each seg, fetch reads, sort them, and create individual count vectors
            for seg in cds_part:
                reads = gnd.get_reads(seg)
                for read in filter(lambda x: len(x.positions) in read_dict,reads):
                    read_dict[len(read.positions)].append(read)
    
                # map and sort by length
                for read_length in read_dict:
                    count_vector = list(gnd.map_fn(read_dict[read_length],seg)[1])
                    count_vectors[read_length].extend(count_vector)
            
            # add each count vector for each length to total
            for k, vec in count_vectors.items():
                counts = numpy.array(vec)
                if cds_part.strand == "-":
                    counts = counts[::-1]
               
                if len(counts) % 3 == 0:
                    counts = counts.reshape((len(counts)/3,3))
                else:
                    if using_roi == False:
                        message = "Length of '%s' coding region (%s nt) is not divisible by 3. Ignoring last partial codon." % (roi.get_name(),len(counts))
                        warnings.warn(message,DataWarning)
                    newlen = len(counts)//3
                    counts = counts[:3*newlen]
                    counts = counts.reshape(newlen,3)
    
                phase_sums[k] += counts[codon_buffer:back_buffer,:].sum(0)

    printer.write("Counted %s ROIs total." % (n+1))
    for k in dtmp:
        dtmp[k] = numpy.array(dtmp[k])
    
    # total reads counted for each size    
    for k in read_lengths:
        dtmp["reads_counted"][dtmp["read_length"] == k] = phase_sums[k].sum() 
    
    # read length distribution
    dtmp["fraction_reads_counted"] = dtmp["reads_counted"].astype(float) / dtmp["reads_counted"].sum() 
    
    # phase vectors
    phase_vectors = { K : V.astype(float)/V.astype(float).sum() for K,V in phase_sums.items() }
    for i in range(3):
        dtmp["phase%s" % i] = numpy.zeros(len(dtmp["read_length"])) 

    for k, vec in phase_vectors.items():
        for i in range(3):
            dtmp["phase%s" % i][dtmp["read_length"] == k] = vec[i]
    
    # phase table
    fn = "%s_phasing.txt" % args.outbase
    printer.write("Saving phasing table to %s ..." % fn)
    dtmp = pd.DataFrame(dtmp)
    with argsopener(fn,args) as fh:
        dtmp.to_csv(fh,columns=["read_length",
                                "reads_counted",
                                "fraction_reads_counted",
                                "phase0",
                                "phase1",
                                "phase2",
                                ],
                       float_format="%.6f",
                       na_rep="nan",
                       sep="\t",
                       index=False,
                       header=True
                      )
        fh.close()
    
    fig = {}
    if args.figsize is not None:
        fig["figsize"] = tuple(args.figsize)

    colors = pp.get_colors_from_args(args,len(read_lengths))

    fn = "%s_phasing.%s" % (args.outbase,args.figformat)
    printer.write("Plotting to %s ..." % fn)
    plot_counts = numpy.vstack([V for (_,V) in sorted(phase_sums.items())])
    fig, (ax1,_) = phase_plot(plot_counts,labels=read_lengths,lighten_by=0.3,
                                cmap=None,color=colors,fig=fig)

    if args.title is not None:
        ax1.set_title(args.title)
    else:
        ax1.set_title("Phasing stats for %s" % args.outbase)

    fig.savefig(fn,dpi=args.dpi,bbox_inches="tight")
    

if __name__ == "__main__":
    main()
