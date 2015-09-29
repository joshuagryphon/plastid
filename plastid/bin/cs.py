#!/usr/bin/env python
"""Count the number of :term:`read alignments<alignment>` and calculate
read densities (in :term:`RPKM`) specifically covering genes.

:term:`Counts <counts>` and densities are calculated separately per gene for
exons, 5' UTRs, coding regions, and 3' UTRs. In addition, positions overlapped
by multiple genes are excluded, as are positions annotated in
:term:`mask annotation files<crossmap>`, if one is provided.

To avoid unnecessarily repeating calculations of gene overlap, the script's
operation is divided into three subprograms:

Generate
    The :func:`generate <do_generate>` mode first takes a gene annotation file,
    finds all genes that share exons, and maps these to a "merged" genes that
    represent their groups. The reason to do this is that in some genome
    annotations, polycistronic messages (such as *polished rice* in
    *D. melanogaster*) are annotated as entirely separate genes, when in reality
    they are a single gene. Positions covered by more than one non-merged gene
    on the same strand are excluded from analysis.

    The remaining positions in each merged gene are then divided
    into the following groups:

        *exon*
            all positions in all transcripts mapping to the merged gene

        *CDS*
            positions which appear in coding regions in *all* transcript
            isoforms mapping to the merged gene. i.e. These positions
            are never part of a fiveprime or threeprime UTR in *any*
            transcript mapping to the merged gene

        *UTR5*
            positions which are annotated only as *5' UTR* in all
            transcript isoforms mapping to the merged gene

        *UTR3*
            positions which are annotated only as *3 UTR* in all
            transcript isoforms mapping to the merged gene

        *masked*
            positions excluded from analyses as directed in an optional
            :term:`mask file`

    .. Rubric :: Output files

    The following files are output, where `OUTBASE` is a name supplied
    by the user:

        OUTBASE_gene.positions
            Tab-delimited text file. Each line is a merged gene, and columns
            indicate the genomic coordinates and lengths of each of the position
            sets above.

        OUTBASE_transcript.positions
            Tab-delimited text file. Each line is a transcript, and columns
            indicate the genomic coordinates and lengths of each of the position
            sets above.

        OUTBASE_gene_REGION.bed
             `BED`_ files showing position sets for `REGION`,
             where `REGION` is one of *exon*, *utr5*, *cds*, *utr3*, or
             *masked*. These contain the same information in
             ``OUTBASE_gene.positions``, but can be visualized easily in a
             :term:`genome browser`


Count
    A :func:`count <do_count>` mode, in which the number of reads mapping to
    each type of position (e.g. *exon*, *CDS*, *UTR5*, and *UTR3*) is tabulated
    for each merged gene. Output is given in a tab-delimited text file that gives
    raw read counts, as well as :term:`RPKM` and corrected lengths for each
    feature.


Chart
    A :func:`chart <do_chart>` mode, which takes output from the ``count`` mode
    and generates several tables and charts that provide broad overviews
    of the data.

See command-line help for each subprogram for details on each mode

See also
--------
:mod:`~plastid.bin.counts_in_region` script
    Count the number of :term:`read alignments <alignment>` covering arbitrary
    regions of interest -- rather than merged genes, as are counted here in
    :mod:`~plastid.bin.cs` -- in the genome, and calculate read densities
    (in reads per nucleotide and in :term:`RPKM`) over these regions.
"""

from plastid.util.scriptlib.argparsers import get_genome_array_from_args,\
                                           get_transcripts_from_args,\
                                           get_alignment_file_parser,\
                                           get_annotation_file_parser,\
                                           get_mask_file_parser,\
                                           get_genome_hash_from_mask_args


from plastid.util.scriptlib.help_formatters import format_module_docstring

from plastid.genomics.roitools import positions_to_segments, SegmentChain
from plastid.genomics.genome_hash import GenomeHash
from plastid.util.io.openers import opener, get_short_name, argsopener
from plastid.util.io.filters import NameDateWriter
from plastid.util.services.sets import merge_sets
from plastid.util.services.decorators import skipdoc
from scipy.misc import comb as combination

import os
import sys
import numpy
import pandas as pd
import scipy.optimize
import argparse
import scipy.stats
import itertools
import copy
import inspect
import gc
import warnings

warnings.simplefilter("once")
printer = NameDateWriter(get_short_name(inspect.stack()[-1][1]))


#===============================================================================
# 'generate' subprogram
#===============================================================================

@skipdoc
def write_output_files(table,title,args):
    """Write gene info table from :py:func:`do_generate` to several output files:

        OUTBASE_gene.positions
            Tab-delimited text file. Each line is a merged gene, and columns
            indicate the genomic coordinates and lengths of each of the position
            sets above.

        OUTBASE_transcript.positions
            Tab-delimited text file. Each line is a transcript, and columns
            indicate the genomic coordinates and lengths of each of the position
            sets above.

        OUTBASE_gene_REGION.bed
             `BED`_ files showing position sets for `REGION`,
             where `REGION` is one of *exon*, *utr5*, *cds*, *utr3*, or
             *masked*. These contain the same information in
             ``OUTBASE_gene.positions``, but can be visualized easily in a
             :term:`genome browser`


    Parameters
    ----------
    table : :class:`pandas.DataFrame`
        Gene info table made in :py:func:`do_generate`

    title : str
        Title ("gene" or "transcript")

    args : :py:class:`argparse.Namespace`
        Command-line arguments
    """
    keys = ("utr5","utr3","cds","masked","exon")
    bed_columns = ["%s_bed" % K for K in keys]

    bedfiles = { X : argsopener("%s_%s_%s.bed" % (args.outbase,title,X),args) for X in keys }

    for _, row in table[bed_columns].iterrows():
        for k in keys:
            bedfiles[k].write(row["%s_bed" % k])

    for k in bedfiles:
        bedfiles[k].close()

    pos_out = argsopener("%s_%s.positions" % (args.outbase,title),args)
    table.to_csv(pos_out,
                 sep="\t",
                 header=True,
                 index=False,
                 na_rep="nan",
                 float_format="%.8f",
                 columns=["region",
                          "exon",
                          "utr5",
                          "cds",
                          "utr3",
                          "masked",
                          "exon_unmasked",
                          "transcript_ids"])
    pos_out.close()

def merge_genes(tx_ivcs):
    """Merge genes whose transcripts share exons into a combined, "merged" gene

    Parameters
    ----------
    tx_ivcs : dict
        Dictionary mapping unique transcript IDs to |Transcripts|


    Returns
    -------
    dict
        Dictionary mapping raw gene names to the names of the merged genes
    """
    dout = {}
    exondicts = {
        "+" : {},
        "-" : {}
    }

    # backmap exons to genes as tuples
    printer.write("Mapping exons to genes...")
    for txid in tx_ivcs.keys():
        my_segmentchain  = tx_ivcs[txid]
        my_gene = my_segmentchain.get_gene()
        chrom  = my_segmentchain.spanning_segment.chrom
        strand = my_segmentchain.spanning_segment.strand
        for my_iv in my_segmentchain:
            start = my_iv.start
            end = my_iv.end
            # separate exons by chromosome and strand to reduce 
            # downstream comparisons in merging
            if chrom not in exondicts[strand]:
                exondicts[strand][chrom] = {}
            try:
                exondicts[strand][chrom][(start,end)].append(my_gene)
            except KeyError:
                exondicts[strand][chrom][(start,end)] = [my_gene]

    for strand in exondicts:
        for chrom in exondicts[strand]:
            exondicts[strand][chrom] = { K : set(V) for K,V in exondicts[strand][chrom].items() }
            printer.write("Flattening genes on %s(%s)..." % (chrom,strand))
            gene_groups = merge_sets(exondicts[strand][chrom].values(),
                                     printer=printer)
 
            for group in gene_groups:
                merged_name = ",".join(sorted(group))
                for gene in group:
                    dout[gene] = merged_name
        
    printer.write("Flattened to %s groups." % len(dout))
    return dout

def process_partial_group(transcripts,mask_hash,printer):
    """Correct boundaries of merged genes, as described in :func:`do_generate`

    Parameters
    ----------
    transcripts : dict
        Dictionary mapping unique transcript IDs to |Transcripts|.
        This set should be complete in the sense that it should contain
        all transcripts that have any chance of mutually overlapping
        each other (e.g. all on same chromosome and strand). 

    mask_hash : |GenomeHash|
        |GenomeHash| of regions to exclude from analysis


    Returns
    -------
    :class:`pandas.DataFrame`
        Table of merged gene positions

    :class:`pandas.DataFrame`
        Table of adjusted transcript positions

    :class:`dict`
        Dictionary mapping raw gene names to merged gene names
    """
    gene_table = { "region"          : [],
                   "transcript_ids"  : [],
                   "exon_unmasked"   : [],
                   "exon"            : [],
                   "masked"          : [],
                   "utr5"            : [],
                   "cds"             : [],
                   "utr3"            : [],
                   "exon_bed"        : [],
                   "utr5_bed"        : [],
                   "cds_bed"         : [],
                   "utr3_bed"        : [],
                   "masked_bed"      : [],
                 }

    # data table for transcripts
    transcript_table = { "region"          : [],
                         "exon"            : [],
                         "utr5"            : [],
                         "cds"             : [],
                         "utr3"            : [],
                         "masked"          : [],
                         "exon_unmasked"   : [],
                         "transcript_ids"  : [],
                         "exon_bed"        : [],
                         "utr5_bed"        : [],
                         "cds_bed"         : [],
                         "utr3_bed"        : [],
                         "masked_bed"      : [],
                       }

    keycombos = list(itertools.permutations(("utr5","cds","utr3"),2))

    # merge genes that share exons & write output
    printer.write("Collapsing genes that share exons...")
    merged_genes = merge_genes(transcripts)

    # remap transcripts to merged genes
    # and vice-versa
    merged_gene_tx = {}
    tx_merged_gene = {}
    printer.write("Mapping transcripts to merged genes...")
    for txid in transcripts:
        my_tx = transcripts[txid]
        my_gene = my_tx.get_gene()
        my_merged = merged_genes[my_gene]
        tx_merged_gene[txid] = my_merged
        try:
            merged_gene_tx[my_merged].append(txid)
        except KeyError:
            merged_gene_tx[my_merged] = [txid]

    # flatten merged genes
    printer.write("Flattening merged genes, masking positions, and labelling subfeatures...")
    for n, (gene_id, my_txids) in enumerate(merged_gene_tx.items()):
        if n % 1000 == 0 and n > 0:
            printer.write("    %s genes..." % n)

        my_gene_positions = []
        chroms  = []
        strands = []
        for my_txid in my_txids:
            my_segmentchain = transcripts[my_txid]
            chroms.append(my_segmentchain.chrom)
            strands.append(my_segmentchain.strand)
            my_gene_positions.extend(my_segmentchain.get_position_list())

            try:
                assert len(set(chroms)) == 1
            except AssertionError:
                printer.write("Skipping gene %s which contains multiple chromosomes: %s" % (gene_id,",".join(chroms)))

            try:
                assert len(set(strands)) == 1
            except AssertionError:
                printer.write("Skipping gene %s which contains multiple strands: %s" % (gene_id,",".join(strands)))

        my_gene_positions = set(my_gene_positions)
        gene_ivc_raw = SegmentChain(*positions_to_segments(chroms[0],strands[0],my_gene_positions))
        gene_table["region"].append(gene_id)
        gene_table["transcript_ids"].append(",".join(sorted(my_txids)))
        gene_table["exon_unmasked"].append(gene_ivc_raw)

    printer.write("    %s genes total." % (n+1))

    # mask genes
    printer.write("Masking positions and labeling subfeature positions...")
    gene_hash = GenomeHash(gene_table["exon_unmasked"],do_copy=False)

    for n,(gene_id,gene_ivc_raw) in enumerate(zip(gene_table["region"],gene_table["exon_unmasked"])):
        if n % 2000 == 0:
            printer.write("    %s genes..." % n)

        my_chrom  = gene_ivc_raw.spanning_segment.chrom
        my_strand = gene_ivc_raw.spanning_segment.strand

        masked_positions = []
        nearby_genes = gene_hash[gene_ivc_raw]

        # don't mask out positions from identical gene
        gene_ivc_raw_positions = gene_ivc_raw.get_position_set()
        nearby_genes = [X for X in nearby_genes if X.get_position_set() != gene_ivc_raw_positions]
        for gene in nearby_genes:
            masked_positions.extend(gene.get_position_list())

        nearby_masks = mask_hash[gene_ivc_raw]
        for mask in nearby_masks:
            masked_positions.extend(mask.get_position_list())

        masked_positions = set(masked_positions)

        gene_positions_raw = gene_ivc_raw.get_position_set()
        mask_ivc_positions = gene_positions_raw & masked_positions
        total_mask_ivc = SegmentChain(*positions_to_segments(my_chrom,my_strand,mask_ivc_positions),
                                      ID=gene_id)
        gene_table["masked"].append(total_mask_ivc)
        gene_table["masked_bed"].append(total_mask_ivc.as_bed())
        
        gene_post_mask = gene_positions_raw - masked_positions
        gene_post_mask_ivc = SegmentChain(*positions_to_segments(my_chrom,my_strand,gene_post_mask),
                                          ID=gene_id)
        gene_table["exon"].append(gene_post_mask_ivc)
        gene_table["exon_bed"].append(gene_post_mask_ivc.as_bed())
    
        masked_positions = total_mask_ivc.get_position_set()
        tmp_positions = { "utr5"  : set(),
                          "cds"   : set(),
                          "utr3"  : set(),
                         }
        txids  = sorted(merged_gene_tx[gene_id])
        chrom  = gene_post_mask_ivc.chrom
        strand = gene_post_mask_ivc.strand

        # pool transcript positions
        for txid in txids:
            transcript = transcripts[txid]
            
            utr5pos = transcript.get_utr5().get_position_set()
            cdspos  = transcript.get_cds().get_position_set()
            utr3pos = transcript.get_utr3().get_position_set()
            
            tmp_positions["utr5"]  |= utr5pos
            tmp_positions["cds"]   |= cdspos
            tmp_positions["utr3"]  |= utr3pos 
        
        # eliminate positions in which CDS & UTRs overlap from each transcript
        for txid in txids:
            transcript = transcripts[txid]
            transcript_positions = {
                          "utr5"  : transcript.get_utr5().get_position_set(),
                          "cds"   : transcript.get_cds().get_position_set(),
                          "utr3"  : transcript.get_utr3().get_position_set(),
                         }

            for key1, key2 in keycombos:
                transcript_positions[key1] -= tmp_positions[key2]
                transcript_positions[key1] -= masked_positions
        
            transcript_table["region"].append(txid)
            
            # all unmasked positions
            my_chain = SegmentChain(*positions_to_segments(chrom,strand,transcript.get_position_set() - masked_positions),
                                    ID=txid)
            transcript_table["exon"].append(str(my_chain))
            transcript_table["exon_bed"].append(my_chain.as_bed())

            # all uniquely-labeled unmasked positions            
            for k,v in transcript_positions.items():
                my_chain = SegmentChain(*positions_to_segments(chrom,strand,v),
                                        ID=txid)
                transcript_table[k].append(str(my_chain))
                transcript_table["%s_bed" % k].append(my_chain.as_bed())
            
            total_mask_ivc.attr["ID"] = txid
            transcript_table["masked"].append(str(total_mask_ivc))
            transcript_table["masked_bed"].append(total_mask_ivc.as_bed())
            transcript_table["exon_unmasked"].append(str(transcript))
            transcript_table["transcript_ids"].append(txid)
        
        tmp_positions2 = copy.deepcopy(tmp_positions)
        for k1,k2 in keycombos:
            tmp_positions[k1] -= tmp_positions2[k2]
            tmp_positions[k1] -= masked_positions

        for k in (tmp_positions.keys()):
            my_chain = SegmentChain(*positions_to_segments(chrom,strand,tmp_positions[k]),
                                    ID=gene_id)
            gene_table[k].append(str(my_chain))
            gene_table["%s_bed" % k].append(my_chain.as_bed())
    
    printer.write("    %s genes total." % (n+1))

    # cast SegmentChains/Transcripts to strings to keep numpy from unpacking them
    conversion_keys = ["exon","utr5","cds","utr3","masked","exon_unmasked"]
    for k in conversion_keys:
        gene_table[k] = [str(X) for X in gene_table[k]]
        transcript_table[k] = [str(X) for X in transcript_table[k]]
    
    gene_df = pd.DataFrame(gene_table)
    gene_df.sort(columns=["region"],inplace=True)

    transcript_df = pd.DataFrame(transcript_table)
    transcript_df.sort(columns=["region"],inplace=True)

    return gene_df, transcript_df, merged_genes

def do_generate(args):
    """Generate gene position files from gene annotations.
    
    1.  Genes whose transcripts share exons are first merged into merged genes 
        using :func:`merge_genes`.
        
    2.  Within merged genes, all positions are classified. All positions are
        included in a set called *exon*. All positions that appear as coding
        regions in all transcripts (i.e. are never part of a 5'UTR or 3'UTR)
        included in a set called *CDS*. Similarly, all positions that appear
        as 5' UTR or 3' UTR in all transcripts are included in sets called
        *UTR5* or *UTR3*, respectively.
    
    3.  Genomic positions that are overlapped by multiple merged genes are
        excluded from the position sets for those genes.
    
    4.  If a :term:`mask file` is supplied, positions annotated in the mask file
        are also excluded
    
    5.  Output is written using :func:`write_output_files`
    
    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        command-line arguments for ``generate`` subprogram
    """
    # variables for transcript <-> merged gene mapping
    transcripts    = {} 
    merged_genes   = {}

    # data table for merged genes
    gene_table = pd.DataFrame({ 
                   "region"          : [],
                   "transcript_ids"  : [],
                   "exon_unmasked"   : [],
                   "exon"            : [],
                   "masked"          : [],
                   "utr5"            : [],
                   "cds"             : [],
                   "utr3"            : [],
                   "exon_bed"        : [],
                   "utr5_bed"        : [],
                   "cds_bed"         : [],
                   "utr3_bed"        : [],
                   "masked_bed"      : [],
                   })

    # data table for transcripts
    transcript_table = pd.DataFrame({
                         "region"          : [],
                         "exon"            : [],
                         "utr5"            : [],
                         "cds"             : [],
                         "utr3"            : [],
                         "exon_bed"        : [],
                         "utr5_bed"        : [],
                         "cds_bed"         : [],
                         "utr3_bed"        : [],
                         "masked"          : [],
                         "exon_unmasked"   : [],
                         "transcript_ids"  : [],
                         "masked_bed"      : [],
                        })

    keycombos = list(itertools.permutations(("utr5","cds","utr3"),2))

    # data
    is_sorted = (args.sorted == True) or \
                (args.tabix == True) or \
                (args.annotation_format == "BigBed")
    source = get_transcripts_from_args(args,printer=printer)
    mask_hash = get_genome_hash_from_mask_args(args)
    
    # loop conditions
    last_chrom = None
    do_loop = True
    z = 0

    # to save memory, we process one chromosome at a time if input file is sorted
    # knowing that at that moment all transcript parts are assembled
    while do_loop == True:
        try:
            tx = next(source)
        except StopIteration:
            do_loop = False

        # if chromosome is completely processed or EOF
        if (is_sorted and tx.spanning_segment.chrom != last_chrom) or do_loop == False:
            if do_loop == True:
                source = itertools.chain([tx],source)

            if last_chrom is not None or do_loop == False:
                printer.write("Merging genes on chromosome/contig '%s'" % last_chrom)
                my_gene_table, my_transcript_table,my_merged_genes = process_partial_group(transcripts,mask_hash,printer)
                gene_table = pd.concat((gene_table,my_gene_table),axis=0)
                transcript_table = pd.concat((transcript_table,my_transcript_table),axis=0)
                merged_genes.update(my_merged_genes)

            del transcripts
            gc.collect()
            del gc.garbage[:]
            transcripts = {}

            # reset last chrom
            last_chrom = tx.spanning_segment.chrom

        # otherwise, remember transcript
        else:
            transcripts[tx.get_name()] = tx
   
    # write output
    printer.write("Writing output...")
    
    merged_fn = "%s_merged.txt" % args.outbase
    number_merged = len(set(merged_genes.values()))
    printer.write("Collapsed %s genes to %s merged groups. Writing to %s" % (len(merged_genes),number_merged,merged_fn))
    fout = argsopener(merged_fn,args,"w")
    for gene,merged_name in sorted(merged_genes.items()):
        fout.write("%s\t%s\n" % (gene,merged_name))
        
    fout.close()
    
    printer.write("Writing gene table and BED files...")
    write_output_files(gene_table,"gene",args)
    printer.write("Writing transcript summary table and BED files...")
    write_output_files(transcript_table,"transcript",args)
    printer.write("Done!")


#===============================================================================
# 'count' subprogram
#===============================================================================

def do_count(args):
    """Count the number of reads falling within a gene (as specified in a
position file made using the `generate` subcommand). Reads are reported
both in raw :term:`counts` and as :term:`RPKM`, and saved to a 
tab-delimited text file.
    
    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        command-line arguments for ``count`` subprogram    
    """
    keys=("exon","utr5","cds","utr3")
    column_order = ["region"]
    with opener(args.position_file) as pos_fh:
        gene_positions = pd.read_table(pos_fh,sep="\t",comment="#",index_col=None,header=0)
        pos_fh.close()

    # read count files
    gnd = get_genome_array_from_args(args,printer=printer,disabled=["normalize"])
    total_counts = gnd.sum()

    printer.write("Dataset has %s counts in it." % total_counts)
    printer.write("Tallying genes...")

    dtmp = { "region" : [] }
    for x in keys:
        for y in ("reads","length","rpkm"):
            label = "%s_%s" % (x, y)
            dtmp[label] = []
            column_order.append(label)
    
    for i,name in enumerate(gene_positions["region"]):
        dtmp["region"].append(name)
        if i % 500 == 0:
            printer.write("Processed %s genes..." % i)
        
        for k in keys:
            ivc = SegmentChain.from_str(gene_positions[k][i])
            total  = sum(ivc.get_counts(gnd))
            length = ivc.length
            rpkm =( 1000 * 1e6 * total / length / total_counts ) if length > 0 else numpy.nan
            dtmp["%s_reads"  % k].append(total)
            dtmp["%s_length" % k].append(length)
            dtmp["%s_rpkm"   % k].append(rpkm)

    fout = argsopener("%s.txt" % args.outbase,args,"w")
    dtmp = pd.DataFrame(dtmp)
    dtmp.to_csv(fout,sep="\t",header=True,index=False,columns=column_order,na_rep="nan",float_format="%.8f")

    fout.close()
    printer.write("Done.")


#===============================================================================
# 'chart' subprogram
#===============================================================================

@skipdoc
def read_count_file(fh,genes_to_include=None):
    """Read rows selected by gene name from a count file
    
    Parameters
    ----------
    fname : str
        Name of output file from cs.py
        
    genes_to_include : list or None
        List of genes to include in output.
        If `None`, all are included
    
    Returns
    -------
    :class:`pandas.DataFrame`

        Count data
    """
    df = pd.read_table(fh,sep="\t",comment="#",index_col=None,header=0)
    import_mask = numpy.array([True if X in genes_to_include else False for X in df["region"]])

    return df[import_mask]

@skipdoc
def get_bin_mask_by_summed_key(rep_a,rep_b,bins,key="exon_reads"):
    """Bins genes into groups based upon the summed counts in samples
    `rep_a` and `rep_b`
    
    Parameters
    ----------
    rep_a : :class:`pandas.DataFrame`
        replicate a, from :func:`read_count_file`
        
    rep_b : :class:`pandas.DataFrame`
        replicate b, from :func:`read_count_file`
        
    bins : :class:`numpy.ndarray`
        sequence of bins, specified as floors
        
    key : str, optional
        quantity on which to perform binning (Default: `exon_reads`)

    Returns
    -------
    dict[bin]
        = numpy boolean mask
    """
    my_sums = rep_a[key] + rep_b[key]
    bins = sorted(bins)
    dout = dict.fromkeys(bins)
    for k in dout.keys():
        dout[k] = numpy.tile(False,len(my_sums))
    for n, my_sum in enumerate(my_sums):
        for i in range(len(bins)-1):
            bin1 = bins[i]
            bin2 = bins[i+1]
            if my_sum >= bin1 and my_sum < bin2:
                dout[bin1][n] = True
        if my_sum >= bins[-1]:
            dout[bins[-1]][n] = True
    return dout

@skipdoc
def get_nonzero_either_mask(vector_a,vector_b):
    """Returns a numpy array of boolean values indicating where values in two
    vectors are both greater than zero.
    
    Parameters
    ----------
    vector_a : numpy.ndarray
        Array of counts or RPKM
        
    vector_b : numpy.ndarray
        Array of counts or RPKM
    
    Returns
    -------
    numpy.ndarray    
        Boolean array that is `True` where both `vector_a` and `vector_b`
        have values greater than zero, and `False` elsewhere.
    """
    return (vector_a > 0) & (vector_b > 0)
                        
@skipdoc
def get_short_samplename(inp):
    """Creates a sample legend label from a sample filename
    by removing everything including and after ``".txt"``

    Parameters
    ----------
    inp : str
        filename
        
    Returns
    -------
    str
    """
    return get_short_name(inp,separator=os.path.sep,terminator=".txt")

@skipdoc
def gaussian(mu,sigma):
    """Returns a function that calculates a Gaussian specified by `mu` and `sigma`
    for an arbitrary X
    
    Parameters
    ----------
    mu : float
        mean of Gaussian
        
    sigma : float
        standard deviation of Gaussian
    
    Returns
    -------
    function
        Gaussian probability density function
    """
    coeff = (2*numpy.pi*(sigma**2))**-0.5
    def f(x):
        return coeff*numpy.e**( -(x-mu)**2 / (2*(sigma**2)) )
    return f

@skipdoc
def fit_gaussian(my_x,my_y):
    """Fits a Gaussian to data
    
    Parameters
    ----------
    my_x : numpy.ndarray
        vector of x-values
        
    my_y : numpy.ndarray
        vector of y-values
        
    Returns
    -------
    float
        mean from Gaussian fit
    
    float
        standard deviation from Gaussian fit
    """
    fitfunc = lambda x, mu, sigma2: (2*numpy.pi*sigma2)**-0.5 * numpy.exp( -(x-mu)**2 / (2*sigma2) )
    p0 = [0.0, 0.05]
    params, _ = scipy.optimize.curve_fit(fitfunc,my_x,my_y,p0=p0[:])
    return tuple(params)
    
def do_chart(args):
    """Produce log-2 fold change histograms and scatter plots for :term:`count`
and :term:`RPKM` values between two samples, as well as a chart plotting
correlation coefficients as a function of summed read counts in both samples 

    Parameters
    ----------
    args : :py:class:`argparse.Namespace`
        command-line arguments for ``chart`` subprogram       
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import numpy.ma as ma
        
    bins = numpy.array(args.bins)    
    scatter_kwargs = { "marker" : "o", "s" : 0.5 }
    
    # read input files
    printer.write("Reading input files: %s..." % ", ".join(args.infiles))
    list_of_genes = set([X.strip() for X in opener(args.list_of_regions)])
    samples = { get_short_samplename(X) : read_count_file(opener(X),list_of_genes)\
                                          for X in args.infiles }
    
    # set up chart colors
    num_colors = combination(len(args.infiles),2)
    if num_colors > 1:
        color_incr = 1.0/(num_colors - 1)
        my_cm = cm.get_cmap("Spectral")
        colors = itertools.cycle([my_cm(X) for X in numpy.arange(0.0,1.0+color_incr,color_incr)])
    else:
        colors = itertools.cycle(["#000000"])
    
    # binkeys = list of keys on which binning will be performed
    binkeys      = tuple(["%s_reads" % k for k in args.regions])
    sample_names = sorted(samples.keys())    
    comparisons  = filter(lambda x: x[0] != x[1],
                          list(itertools.combinations(sample_names,2)))
    
    comparison_labels      = sorted(["%s_vs_%s" % (X,Y) for X,Y in comparisons]) 
    bigtable               = {}.fromkeys(comparison_labels)
    corrcoef_by_bin_table  = {}.fromkeys(comparison_labels)
    
    # ki, kj = names of samples i and j
    # vi, vj = data of samples i and j    
    for ki, kj in comparisons:
        try:
            assert (samples[ki]["region"] == samples[kj]["region"]).all()
        except AssertionError:
            printer.write("Mismatched line entries for samples %s and %s." % (ki,kj))
        vi = samples[ki]
        vj = samples[kj]
        printer.write("Comparing %s to %s:" % (ki,kj))
        
        label = "%s_vs_%s" % (ki,kj)
        
        bigtable[label]              = {}
        corrcoef_by_bin_table[label] = {}

        for binkey in binkeys:
            bigtable[label][binkey]              = { K : copy.deepcopy({}) for K in args.metrics }
            corrcoef_by_bin_table[label][binkey] = { K : copy.deepcopy({}) for K in args.metrics }
        
        for region in args.regions:
            region_counts = "%s_reads" % region            
            summed_counts = vi[region_counts] + vj[region_counts]
            
            for metric in args.metrics:                                        
                region_metric = "%s_%s" % (region,metric)
                printer.write("    -%s across %s for all >=128..." % (metric,region))
                
                # divide into regions >=128 and <128 an calcu
                vi_128plus = vi[region_metric][summed_counts >= 128]
                vj_128plus = vj[region_metric][summed_counts >= 128]
                pearsonr  = scipy.stats.pearsonr(vi_128plus,vj_128plus)[0]
                spearmanr = scipy.stats.spearmanr(vi_128plus,vj_128plus)[0]
                
                # need to mask out zeros before plotting, because matplotlib
                # assigns zeros values in order to (erroneously) plot on loglog
                # also need to remove zeros for log2 fold-changes (below)

                nonzero_mask = get_nonzero_either_mask(vi[region_metric],vj[region_metric])
                
                vi_128plus_nonzero  = vi[region_metric][(summed_counts >= 128) & nonzero_mask]
                vj_128plus_nonzero  = vj[region_metric][(summed_counts >= 128) & nonzero_mask]
                vi_128minus_nonzero = vi[region_metric][(summed_counts < 128) & nonzero_mask]
                vj_128minus_nonzero = vj[region_metric][(summed_counts < 128) & nonzero_mask]

                # scatter plot
                scatter_title = "%s vs %s (%s %s)" % (ki,kj,region,metric)
                plt.figure()                    
                plt.loglog()
                plt.scatter(vi_128minus_nonzero,vj_128minus_nonzero,
                            color="#FF0000",label="< 128 reads",**scatter_kwargs)                    
                plt.scatter(vi_128plus_nonzero,vj_128plus_nonzero,
                            color="#000000",label=">= 128 reads",**scatter_kwargs)                    

                text_kwargs = { "horizontalalignment" : "left",
                                "transform" : plt.gca().transAxes }
                plt.text(0.02,0.9,"spearman r**2 >= 128: %0.05f" % spearmanr**2,**text_kwargs)
                plt.text(0.02,0.86,"pearson rho**2 >= 128: %0.05f" % pearsonr**2,**text_kwargs)
                
                plt.xlabel(ki)
                plt.ylabel(kj)
                plt.title(scatter_title)                    
                plt.legend(loc="lower right")
                plt.savefig("%s_scatter_%s_%s_%s.%s" % (args.outbase, label, region, metric, args.format))
                plt.close()
                
                
                # log2 fold-change histogram
                log2_ratios = numpy.log2(vj_128plus_nonzero/vi_128plus_nonzero)
                log2_mean   = numpy.mean(log2_ratios)
                log2_std    = numpy.std(log2_ratios)
                min_diff    = 2**(3*log2_std)
                
                plt.figure()
                n,hbins,_ = plt.hist(log2_ratios,bins=50,normed=True)
                hbins_centered = numpy.convolve(hbins,numpy.array([0.5,0.5]),"valid")
                
                # fit a gaussian
                fit_mu,fit_sigma2 = fit_gaussian(hbins_centered,n)
                fit_min_diff = 2**(3*fit_sigma2**0.5)
                plt.bar(hbins[:-1],n,width=hbins[1]-hbins[0],align='edge')
                my_min = numpy.min(log2_ratios)
                my_max = numpy.max(log2_ratios)
                
                # gaussian from calculated parameters                
                x =  numpy.arange(my_min,my_max,0.01)
                norm = [gaussian(log2_mean,log2_std)(X) for X in x]
#                normfit = [gaussian(fit_mu,fit_sigma2**0.5)(X) for X in x]

                # plot both
                plt.plot(x,norm,'r-',linewidth=1,label="Sample")
#                plt.plot(x,normfit,'r.',linewidth=1,label="Fit")
                plt.legend()
                
                text_kwargs = { "horizontalalignment" : "left",
                                "transform" : plt.gca().transAxes }
                plt.text(0.02,0.9, "sample mean : %0.5f" % log2_mean,**text_kwargs)
                plt.text(0.02,0.87,"sample stdev : %0.5f" % log2_std,**text_kwargs)
                plt.text(0.02,0.81,"3-std fold change: %0.2f" % min_diff,**text_kwargs)
#                plt.text(0.02,0.78, "fit mean : %0.5f" % fit_mu,**text_kwargs)
#                plt.text(0.02,0.75,"fit stdev : %0.5f" % fit_sigma2**0.5,**text_kwargs)
#                plt.text(0.02,0.72,"fit fold change : %0.5f" % fit_min_diff,**text_kwargs)                    
                plt.text(0.02,0.69,"genes counted : %s" % len(log2_ratios),**text_kwargs)
                plt.title("log2-fold %s %s change in %s/%s" % (region, metric, kj, ki))
                plt.xlabel("log2-fold change of %s/%s" % (kj,ki))
                plt.ylabel("Frequency")
                plt.savefig("%s_log2hist_%s_%s_%s.%s" % (args.outbase, label, region, metric, args.format))
                
                
                # add entries to bigtable for export later
                bigtable[label][region_counts][metric]["pearsonr"]          = pearsonr
                bigtable[label][region_counts][metric]["spearmanr"]         = spearmanr
                bigtable[label][region_counts][metric]["log2_mean"]         = log2_mean
                bigtable[label][region_counts][metric]["log2_std"]          = log2_std
                bigtable[label][region_counts][metric]["min_diff"]          = min_diff
                bigtable[label][region_counts][metric]["fit_log2_mean"]     = fit_mu
                bigtable[label][region_counts][metric]["fit_log2_std"]      = fit_sigma2**0.5
                bigtable[label][region_counts][metric]["fit_min_diff"]      = fit_min_diff
                bigtable[label][region_counts][metric]["num_genes_128plus"] = len(vi_128plus)
                bigtable[label][region_counts][metric]["num_genes_log2"]    = len(vi_128plus_nonzero)

                
                # do bin-by-bin counting
                printer.write("    -%s across %s by bin..." % (metric,region))                    
                bin_masks = get_bin_mask_by_summed_key(vi,vj,bins=args.bins,key=region_counts)
                for my_bin,bin_mask in sorted(bin_masks.items()):
                    
                    bin_vec_i = vi[region_metric][bin_mask]
                    bin_vec_j = vj[region_metric][bin_mask]
                    
                    # make sure there are genes in bin before attempting calculations
                    if len(bin_vec_i) > 0:
                        nonzero_binmask = get_nonzero_either_mask(bin_vec_i,bin_vec_j)
                        bin_vi_log2 = bin_vec_i[nonzero_binmask]
                        bin_vj_log2 = bin_vec_j[nonzero_binmask]
                        my_logs      = numpy.log2(bin_vj_log2/bin_vi_log2)
                        my_logmean   = numpy.mean(my_logs)
                        my_logstd    = numpy.std(my_logs)
                        if len(bin_vec_i) > 2:
                            my_pearsonr  = scipy.stats.pearsonr(bin_vec_i,bin_vec_j)[0]
                            my_spearmanr = scipy.stats.spearmanr(bin_vec_i,bin_vec_j)[0]
                        else:
                            my_spearmanr = numpy.nan
                            my_pearsonr = numpy.nan
                    else:
                        # fill with dummy values
                        my_logs      = numpy.array([])
                        my_logmean   = numpy.nan
                        my_logstd    = numpy.nan
                        my_spearmanr = numpy.nan
                        my_pearsonr  = numpy.nan
                    
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin] = {}
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["pearsonr"]     = my_pearsonr
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["spearmanr"]    = my_spearmanr
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["genes_in_bin"] = sum(bin_mask)
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["log2_genes_in_bin"] = sum(nonzero_binmask)
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["log2_mean"]    = my_logmean
                    corrcoef_by_bin_table[label][region_counts][metric][my_bin]["log2_std"]     = my_logstd

    # export big (non-binned) table
    printer.write("Writing tables...")
    
    bigtable_out = argsopener("%s_bigtable.txt" % args.outbase,args,"w")
    stats = ("num_genes_128plus",  # 0
             "pearsonr",           # 1
             "spearmanr",          # 2
             "num_genes_log2",     # 3
             "log2_mean",          # 4
             "log2_std",           # 5
             "min_diff",           # 6
             "fit_log2_mean",      # 7 
             "fit_log2_std",       # 8
             "fit_min_diff")       # 9
    
    header = ["#region","metric","statistic"]
    header += [X for X in comparison_labels]
    
    bigtable_out.write("\t".join(header)+"\n")
    for region in binkeys:
        for metric in args.metrics:
            for stat in stats:
                ltmp = [region,metric,stat]
                for key in comparison_labels:
                    ltmp.append(bigtable[key][region][metric][stat])
                ltmp = [str(X) for X in ltmp]
                bigtable_out.write("\t".join(ltmp)+"\n")
    
    bigtable_out.close()
    
    
    # export binned data table and make bin-by-bin plots
    bintable_out = argsopener("%s_bintable.txt" % args.outbase,args,"w")
    
    region_metrics = ["%s_%s" % (X,Y) for X in args.regions for Y in args.metrics]
    stats = ["genes_in_bin",           # 0
             "pearsonr",               # 1
             "spearmanr",              # 2
             "",                       # 3
             "log2_genes_in_bin",      # 4
             "log2_mean",              # 5
             "log2_std",               # 6
             ""]                       # 7
    for region_metric in region_metrics:
        plt.close()
        plt.figure()
        plt.semilogx(basex=2)
        plt.title("Correlation coefficients by bin for %s" % region_metric)
        plt.xlabel("Bin")
        plt.ylabel("Spearman rho")

        region,metric = region_metric.split("_")
        region_counts = "%s_reads" % region
        bintable_out.write("%s\n\n" % region_metric)
        
        for label in comparison_labels:
            corrcoefs = []
            bintable_out.write("%s\t\t\t\t\t\t\t\t" % label)
            bintable_out.write("\n")
            bintable_out.write("bin\t" + "\t".join(stats)+"\n")
            for my_bin in bins:
                ltmp = [my_bin]
                for stat in stats:
                    ltmp.append(corrcoef_by_bin_table[label][region_counts][metric][my_bin].get(stat,""))
                ltmp = [str(X) for X in ltmp]
                bintable_out.write("\t".join(ltmp)+"\n\n")
                corrcoefs.append(corrcoef_by_bin_table[label][region_counts][metric][my_bin]["spearmanr"])
            corrcoefs = ma.masked_invalid(corrcoefs)
            plt.plot(bins[~corrcoefs.mask],
                     corrcoefs[~corrcoefs.mask],
                     label=label,
                     color=next(colors))
        plt.legend(loc="lower right")
        plt.savefig("%s_corrcoef_by_bin_%s_%s.%s" % (args.outbase,region_metric,label,args.format))

        bintable_out.write("\n\n")
    bintable_out.close()
    
    printer.write("Done.")
        


#===============================================================================
# main program body
#===============================================================================

def main(argv=sys.argv[1:]):
    """Command-line program
    
    Parameters
    ----------
    argv : list, optional
        A list of command-line arguments, which will be processed
        as if the script were called from the command line if
        :func:`main` is called directly.

        Default: `sys.argv[1:]`. The command-line arguments, if the script is
        invoked from the command line
    """

    alignment_file_parser  = get_alignment_file_parser(disabled=["normalize"])
    annotation_file_parser = get_annotation_file_parser()
    mask_file_parser = get_mask_file_parser()
    
    generator_help = "Create unambiguous position file from GFF3 annotation"
    generator_desc = format_module_docstring(do_generate.__doc__)

    counter_help   = "Count reads in unambiguous gene positions"
    counter_desc   = format_module_docstring(do_count.__doc__)

    chart_help     = "Produce charts comparing reads between samples"
    chart_desc     = format_module_docstring(do_chart.__doc__)

    parser = argparse.ArgumentParser(description=format_module_docstring(__doc__),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title="subcommands",
                                       description="choose one of the following",
                                       dest="program")
    gparser    = subparsers.add_parser("generate",
                                       help=generator_help,
                                       description=generator_desc,
                                       formatter_class=argparse.RawDescriptionHelpFormatter,
                                       parents=[annotation_file_parser,
                                                mask_file_parser],
                                       )
    cparser    = subparsers.add_parser("count",
                                       help=counter_help,
                                       description=counter_desc,
                                       parents=[alignment_file_parser],
                                       formatter_class=argparse.RawDescriptionHelpFormatter,
                                       )
    pparser    = subparsers.add_parser("chart",
                                       help=chart_help,
                                       description=chart_desc,
                                       formatter_class=argparse.RawDescriptionHelpFormatter)

    gparser.add_argument("outbase",metavar="outbase",type=str,
                         help="Basename for output files")
    
    cparser.add_argument("position_file",type=str,metavar="file.positions",
                         help="File assigning positions to genes or transcripts (made using 'generate' subcommand)")
    cparser.add_argument("outbase",type=str,help="Basename for output files")
    
    pparser.add_argument("-i","--in",nargs="+",type=str,dest="infiles",
                         help="input files, made by 'count' subprogram")
    pparser.add_argument("--bins",nargs="+",type=int,
                         default=(0,32,64,128,256,512,1024,2048,4096),
                         help="Bins into which features are partitioned based on counts")
    pparser.add_argument("--regions",nargs="+",type=str,
                         default=("exon","utr5","cds","utr3"),
                         help="Regions to compare (default: exon, utr5, cds, utr3)")
    pparser.add_argument("--metrics",nargs="+",type=str,
                         default=("rpkm","reads"),
                         help="Metrics to compare (default: rpkm, reads)")
    pparser.add_argument("--format",default="pdf",choices=("pdf","svg","png"),
                         help="Output format for charts (default: pdf)")
    pparser.add_argument("list_of_regions",type=str,metavar='gene_list.txt',
                         help="File listing regions (genes or transcripts), one per line, to include in count")
    pparser.add_argument("outbase",type=str,
                         help="Basename for output files")

    args = parser.parse_args(argv)

    if args.program == "generate":
        #generate position file
        do_generate(args)
    
    elif args.program == "count":
        #use position file to count gene expression in infiles
        do_count(args)

    elif args.program == "chart":
        #use count files to generate a family of charts and tables
        do_chart(args)
        
if __name__ == "__main__":
    main()
