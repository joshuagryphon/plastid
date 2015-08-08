Performing metagene analyses
============================

 .. TODO update graphics to reflect 5'/14 mapping

Definition
----------
A :term:`metagene` analysis is an average of quantitative data over one or more
genomic regions (often genes or transcripts) aligned at some internal feature.
For example, a :term:`metagene` profile could be built around:

  - the average of ribosome density surrounding the start codons of all 
    transcripts in a :term:`ribosome profiling` dataset
  
  - an average phylogenetic conservation score surounding the 5' splice site of
    the first introns of all transcripts

Often, this sort of analysis reveals patterns across regions that may not be
obvious when looking at any individual region.


Scope of this tutorial
----------------------
In this tutorial, we:

  - :ref:`Use the metagene script <metagene-script>` to perform a :term:`metagene`
    analysis of ribosome density surrounding start codons in a
    :term:`ribosome profiling` dataset

  - Discuss how to write functions to 
    :ref:`perform custom metagene analyses surrounding arbitary features <concepts-metagene-roll-your-own>`
    of interest. As an example, we make a metagene average of ribosome density
    :ref:`surrounding the highest peak <metagene-big-spike>` of ribosome density
    in coding regions.

This tutorial uses the :doc:`/test_dataset` we have prepared.


Computation
-----------

A metagene average is computed by:

 #. Fetching vectors of *quantitative data* -- in this example,
    :term:`ribosome-protected footprints` -- surrounding each occurrence
    of the *feature of interset* -- in this case, a start codon -- in each
    region:

     .. figure:: /_static/images/metagene_unnormalized_vectors.png
        :alt: Unnormalized vectors of ribosome counts
        :figclass: captionfigure
        
        Unnormalized vectors of :term:`ribosome-protected footprints <footprint>`,
        shown above corresponding transcript models. Thick boxes: coding
        regions. Thin boxes: 5' UTRs.


 #. Normalizing each vector to the same scale by dividing by the total number
    of counts in a normalization window (here, the last 30 nucleotides of
    the portion of coding region shown):

     .. image:: /_static/images/metagene_normalized_vectors.png
        :alt: Normalized vectors of ribosome counts

 
 #. Aligning each vector at the feature of interest (the start codon)
 
 #. Taking an average (e.g. a median, mean, or percentile) over all of the
    normalized vectors at each aligned nucleotide position. Typically,
    we use the median:

     .. figure:: /_static/images/metagene_average_profile.png
        :alt: Creation of metagene profile
        :figclass: captionfigure

        Final :term:`metagene` average over the four transcripts shown above.
             

 .. _metagene-script:
 
Using the |metagene| script
---------------------------
When called from the command-line, |metagene| divides the analysis into 
two subprograms:

 - :ref:`metagene-generate` processes a genome annotation into a set 
   of aligned windows surrounding a feature of interest (e.g. a start codon)
 
 - :ref:`metagene-count` takes the revised annotation from :ref:`generate <metagene-generate>`
   and a dataset of :term:`read alignments` or other quantitative data,
   and produces the metagene average
  
For convenience, a :ref:`chart <metagene-chart>` subprogram is also provided to
plot the outcome of one or more runs of the :ref:`count <metagene-count>`
subprogram.
 
 .. _metagene-generate:

The ``generate`` subprogram
...........................
The first step in a |metagene| analysis is to examine all of the regions of
interest -- in this example, transcripts -- in a genome annotation, detect an
interesting sub-feature -- here, a start codon -- and build a sub-window within
the transcript that surrounds the sub-feature. To do so, the |metagene|
:ref:`generate <metagene-generate>` subprogram performs the following steps:

 #. Group transcripts by gene, so each start codon is only counted once.

 #. If a gene has multiple start codons, exclude that gene.
    
    If a gene has a single start codon, then find the
    :term:`maximal spanning window` of the gene, defined as the largest possible
    window surrounding the start codon in which **all** transcripts from
    that gene map to **identical** genomic positions. This prevents any
    positional ambiguity from entering the average:

    .. figure:: /_static/images/metagene_maximal_spanning_window.png
       :alt: Metagene - maximal spanning window
       :figclass: captionfigure

 #. Once a :term:`maximal spanning window` is defined for a gene, determine the location
    of the start codon relative to the window, so that the maximal spanning
    windows for all genes may be aligned at the start codon during the 
    :ref:`count <metagene-count>` step.
 
 #. Save the maximal spanning windows to a `BED`_ file for inspection in a
    :term:`genome browser` or other analysis pipeline, and to a
    text file (called an *ROI file*) for use in the :ref:`count <metagene-count>`
    subprogram.

The :ref:`generate <metagene-generate>` program only needs to be run once per landmark of interest
(i.e. once for stop codons, once for start codons), and only needs to be 
re-run if the genome :term:`annotation` changes (e.g. due to revisions,
additions, or deletions of gene/transcript models).

We will run the :ref:`generate <metagene-generate>` program using the following arguments:

  - ``merlin_cds_start`` : name all output files with the prefix `merlin_cds_start` 
  - ``--landmark cds_start`` : calculate a metagene average surrounding start codons
  - ``--annotation_files merlin_orfs.gtf`` : use transcript models from this annotation file
  - ``--downstream 200`` : include up to 200 nucleotides downstream of the start codon in the average

The program is called from the terminal: 

 .. code-block:: shell

    $ metagene generate merlin_cds_start --landmark cds_start --annotation_files merlin_orfs.gtf --downstream 200

For a detailed description of these and other command-line arguments, see the
:mod:`metagene script documentation <yeti.bin.metagene>`


 .. _metagene-count:
 
The ``count`` subprogram
........................
Once :ref:`generate <metagene-generate>` has made an ROI file, |metagene|'s
:ref:`count <metagene-count>` subprogram can be used to tabulate metagene averages.
Specifically, :ref:`count <metagene-count>` performs the following steps:

 #. For each :term:`maximal spanning window` in the ROI file:

     #. fetch a vector of counts at each position from the sample dataset
        (in this case, :term:`ribosome profiling` alignments).

     #. If the vector of counts has a sufficient number of alignments within
        a user-specified *normalization region*, include it. Otherwise, exclude
        the vector.
     
     #. Normalize the vector by the number of counts in the *normalization region.*

 #. Construct a metagene average by taking the median over all normalized vectors
    at each position, excluding any vectors that happen not cover that position
    (e.g. because the maximal spanning window for that gene was too small).

 #. For each position, save the metagene average and the number of genes included
    in the average to a tab-delimited text file.

To call the :ref:`count <metagene-count>` program, type into a terminal window.
In this example ``--fiveprime --offset 14`` specify our :term:`mapping rule` for the
:term:`P-site offset`, estimating the offset as 14 nucleotides from the 5' end
of each :term:`read alignment`.

 .. code-block:: shell

    $ metagene count merlin_cds_start_rois.txt SRR609197_riboprofile --count_files SRR609197_riboprofile.bam --fiveprime --offset 14


A number of files are created and may be used for further processing. In our
example, ``SRR609197_riboprofile_metagene_profile.txt`` contains the final, reduced data.
This file contains three columns:

  #. *x:* an X-coordinate indicating the distance in nucleotides from
     the start codon
  
  #. *metagene_profile:* the value of the metagene average at *x*
       
  #. *regions_counted:* the number of regions included in the average at *x*


 .. _metagene-chart:

The ``chart`` subprogram
........................
For convenience, a :ref:`chart <metagene-chart>` subprogram is included. It can plot multiple
metagene profiles (each from a run of the :ref:`count <metagene-count>` subprogram) in
a single plot:

 .. code-block:: shell
 
    $ metagene chart SRR609197_riboprofile_cds_start.png  SRR609197_riboprofile_metagene_profile.txt --landmark "start codon" --title "Metagene demo"

This produces the image:

 .. figure:: /_static/images/demo_metagene_cds_start.png
    :align: center
    :alt: metagene profile surrounding start codon
    :figclass: captionfigure
    
    :term:`Metagene` profile surrounding start codon.


.. _concepts-metagene-roll-your-own:

Beyond start and stop codons: defining your own window functions
----------------------------------------------------------------

:term:`Metagene averages <metagene average>` can be useful for other questions,
types of regions, and experimental data. For this reason, |metagene| offers tools
to create maximal spanning windows surrounding any feature of interest.

Window functions
................
To make :term:`maximal spanning windows <maximal spanning window>` around a
feature, |metagene| requires a
*window function*. The *window function* must identify and build a window around
the feature of interest (e.g. a start codon) in each individual region examined
(for example, each transcript).

Regions that define a `"gene_id"` attribute in their `attr` dictionaries,
will be grouped by shared `"gene_id"` to make a single
:term:`maximal spanning window` for each group.

|metagene| comes with two window functions:

  - :func:`~yeti.bin.metagene.window_cds_start`, for defining windows
    surrounding start codons

  - :func:`~yeti.bin.metagene.window_cds_stop`, for defining windows
    surrounding stop codons

Once you have defined a window function, :func:`yeti.bin.metagene.do_generate`
can use it to generate :term:`maximal spanning windows <maximal spanning window>`.

Parameters
..........
Window functions must take the following parameters, in order:

    `roi` : |SegmentChain|
        Input ROI for which a window will be generated.
        If `"gene_id"` is defined in `roi.attr`, then
        all `roi`s sharing the same `"gene_id"` will
        be used to generate a single
        :term:`maximal spanning window` covering all of them.

    `flank_upstream` : ``int``
        Nucleotide length upstream of the feature of interest
        to include in the :term:`maximal spanning window`, if `roi` has
        such a feature

    `flank_downstream` : ``int``
        Nucleotide length downstream of the feature of interest
        to include in the :term:`maximal spanning window`, if `roi` has
        such a feature
    
    `ref_delta` : ``int``, optional
        Offset in nucleotides from the *feature of interest* to 
        the *reference point* at which all :term:`maximal spanning window`
        count vectors will be aligned when the metagene average
        is calculated. If `0`, the feature of interest is the
        reference point. (Default: `0`)

Return values
.............
Window functions must return the following values, in order:

    |SegmentChain|
        Window surrounding feature of interest if `roi` has such a feature.
        Otherwise, return a zero-length |SegmentChain|. :func:`yeti.bin.metagene.do_generate`
        will use these to generate maximal spanning windows.
        
    
    ``int``
        offset to align window with all other windows, calculated as
        :samp:`max(0,{flank_upstream} - {distance to reference point from 5' end of window})`,
        or :obj:`numpy.nan` if `roi` does not contain a feature of interest 

    (``str``, ``int``, ``str``)
        Genomic coordinate of reference point as `(chromosome name, coordinate, strand)`
        or :obj:`numpy.nan` if `roi` does not contain a feature of interest

Here is a window function that produces windows surrounding transcription
start sites::

    >>> def window_transcript_start(roi,flank_upstream,flank_downstream,ref_delta=0):
    >>>     """Window function for metagenes surrounding transcription start sites
    >>> 
    >>>     Returns
    >>>     -------
    >>>     SegmentChain
    >>>         Window surrounding transcript start site
    >>> 
    >>>     int
    >>>         Offset to align window with all other windows
    >>> 
    >>>     (str,int,str)
    >>>         Genomic coordinate of transcription start site as *(chromosome name, coordinate, strand)*
    >>>     """
    >>>     chrom,tx_start_genome,strand = roi.get_genomic_coordinate(0)
    >>>     segs = roi.get_subchain(0,flank_downstream)
    >>> 
    >>>     if strand == "+":
    >>>         new_segment_start = tx_start_genome - flank_upstream
    >>>         # need to add one for half-open coordinate end positions
    >>>         new_segment_end = roi.get_genomic_coordinate(flank_downstream)[1] + 1
    >>>         offset = 0
    >>>     else:
    >>>         new_segment_start = roi.get_genomic_coordinate(flank_downstream)[1]
    >>>         # need to add one for half-open coordinate end positions
    >>>         new_segment_end = tx_start_genome + flank_upstream + 1
    >>>         if roi.get_length() < flank_downstream:
    >>>             offset = flank_downstream - roi.get_length()
    >>>         else:
    >>>             offset = 0
    >>> 
    >>>     outside_segment = GenomicSegment(chrom,
    >>>                                      new_segment_start,
    >>>                                      new_segment_end,
    >>>                                      strand)
    >>>     segs.add_segments(outside_segment)
    >>>     new_chain = SegmentChain(*tuple(segs))
    >>> 
    >>>     # ref point is always `flank_upstream` from window in this case
    >>>     return new_chain, offset, new_chain.get_genomic_coordinate(flank_upstream)


Here is a window function that produces windows surrounding the highest spike
in read density in a transcript. Note, it uses data structures in the global
scope:

 .. code-block:: python
 
    >>> import numpy

    >>> def window_biggest_spike(roi,flank_upstream,flank_downstream,ref_delta=0):
    >>>     """Window function for metagenes surrounding peaks of read density
    >>>     
    >>>     ALIGNMENTS must be defined in global scope as a GenomeArray
    >>>     
    >>>     
    >>>     Returns
    >>>     -------
    >>>     SegmentChain
    >>>         Window surrounding transcript start site
    >>> 
    >>>     int
    >>>         Offset to align window with all other windows
    >>> 
    >>>     (str,int,str)
    >>>         Genomic coordinate of transcription start site as *(chromosome name, coordinate, strand)*
    >>>     """
    >>>     counts      = roi.get_counts(ALIGNMENTS)
    >>>     if len(counts) > 0:
    >>>         # ignore first 5 and last 5 codons, which will have big
    >>>         # initiation/termination peaks that we want to skip
    >>>         max_val_pos = 15 + counts[15:-15].argmax()
    >>>         ref_point_genome = roi.get_genomic_coordinate(max_val_pos)
    >>> 
    >>>         new_chain_start = max_val_pos - flank_upstream
    >>>         new_chain_end   = max_val_pos + flank_downstream
    >>> 
    >>>         offset    = 0
    >>> 
    >>>         # if new start is outside region, set it to zero and memorize offset
    >>>         if new_chain_start < 0:
    >>>             offset = -new_chain_start
    >>>             new_chain_start = 0
    >>> 
    >>>         # if new end is outside region, set it to end of region
    >>>         if new_chain_end > roi.get_length():
    >>>             new_chain_end = roi.get_length()
    >>> 
    >>>         new_chain = roi.get_subchain(new_chain_start,new_chain_end)
    >>>     else:
    >>>         new_chain = SegmentChain()
    >>>         offset = ref_point_genome = numpy.nan
    >>>         
    >>>     return new_chain, offset, ref_point_genome

 .. _metagene-big-spike:
 
Once your window function is written, you can generate maximal spanning windows.
Here we use the ``window_biggest_spike()`` function we just wrote::

    >>> from yeti.bin.metagene import do_generate
    >>> from yeti.genomics.genome_hash import GenomeHash
    >>> from yeti.readers.gff import GTF2_TranscriptAssembler

    >>> import pysam
    >>> import numpy
    >>> from yeti.genomics.genome_array import BAMGenomeArray, ThreePrimeMapFactory

    >>> # window_biggest_spike() needs read alignments stored in a variable
    >>> # called ALIGNMENTS. so let's load some
    >>> ALIGNMENTS = BAMGenomeArray([pysam.Samfile("SRR609197_riboprofile.bam","rb")])
    >>> ALIGNMENTS.set_mapping(VariableFivePrimeMapFactory(offset_dict))

    >>> # skip masking out any repetitive regions for purpose of demo
    >>> dummy_mask_hash = GenomeHash([])

    >>> #load features, in our case, transcripts
    >>> transcripts = list(GTF2_TranscriptAssembler(open("merlin_orfs.gtf")))

    >>> # include 100 nucleotides up- and downstream of feature
    >>> flank_upstream = flank_downstream = 100

    >>> data_table, segment_chains = do_generate(transcripts,dummy_mask_hash,
    >>>                                          flank_upstream,flank_downstream,
    >>>                                          window_func=window_biggest_spike)


:meth:`~yeti.bin.metagene.do_generate` returns an |ArrayTable| (similar to a :class:`pandas.DataFrame`)
of data and a list of |SegmentChains| corresponding to the
:term:`maximal spanning windows <maximal spanning window>` for each row
in the |ArrayTable|. It is useful to save these in the same formats that the 
:ref:`generate <metagene-generate>` program uses::

    >>> with open("SRR609197_riboprofile_big_spike_roi_file.txt","w") as roi_fh:
    >>>     data_table.to_file(roi_fh)
    >>>     roi_fh.close()
    >>>     
    >>> with open("SRR609197_riboprofile_big_spike_rois.bed","w") as bed_fh:
    >>>     for roi in segment_chains:
    >>>         bed_fh.write(roi.as_bed())
    >>> 
    >>>     bed_fh.close()
 
Then, you can use the ROI file you just made with |metagene|
:ref:`count <metagene-count>` and :ref:`chart <metagene-chart>` as previously:

 .. code-block:: shell

     $ metagene count SRR609197_riboprofile_big_spike_roi_file.txt SRR609197_riboprofile_big_spike --count_files SRR609197_riboprofile.bam --fiveprime_variable --offset SRR609197_riboprofile_p_offsets_adjusted.txt
     
     $ metagene chart SRR609197_riboprofile_big_spike.png SRR609197_riboprofile_big_spike_metagene_profile.txt --landmark "highest ribosome peak" --title "Custom metagene demo"

Which yields:

 .. figure:: /_static/images/metagene_big_spike_demo.png
    :align: center
    :alt: metagene profile surrounding biggest peak
    :figclass: captionfigure
    
    :term:`Metagene` profile surrounding largest spike of ribosome density in coding
    region, excluding :term:`start <start codon peak>` and
    :term:`stop codon peaks <stop codon peak>`.
    
    

See also
--------
  - Module documentation for |metagene| program for detailed description
    of command-line arguments and output files of the three subprograms

    
