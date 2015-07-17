Performing metagene analyses
============================

Definition
----------
A :term:`metagene` analysis is a position-wise average of quantitative data
over one or more genomic loci (often genes or transcripts). The purpose
is to determine whether, when data is pooled, any interesting pattern emerges.
This is similar to the concept of a spike-triggered average in neuroscience.

For example, in a :term:`ribosome profiling` data set, one might fetch
vectors of ribosome densities at each position of each transcript,
and align these vectors at their start codons. A :term:`metagene average`
is computed by normalizing these vectors to the same scale, and then
taking an average value (e.g. a median, mean, or percentile) of the
aligned, normalized vectors at each position:

    [TODO: graphics]

.. TODO : make graphical image. vectors -> alignment -> normalization -> average


.. TODO : poach graphic from Ingolia et al


Scope of this tutorial
----------------------
In this tutorial, we:

  - Use the |metagene| script to perform a :term:`metagene` analysis
    of ribosome density surrounding start codons in a :term:`ribosome profiling`
    dataset

  - Discuss with examples how to write functions to perform custom metagene
    analyses surrounding arbitary features of interest in a genome or
    transcriptome, for example:

      - A nucleotide sequence occuring in a transcript

      - A peak in read density from a :term:`high-throughput sequencing`
        experiment (:term:`RNA-Seq`, :term:`ClIP-Seq`, :term:`DMS-Seq`)


This tutorial uses the :doc:`/test_dataset` we have prepared.


Using the |metagene| script
---------------------------
|metagene| automates metagene analysis. By default, it performs analyses
surrounding start and stop codons of coding transcripts. However, the 
script is easily extended to perform arbitrary analyses (described
:ref:`below <concepts-metagene-roll-your-own>`). Here, we will stick
to start codons.

The ``generate`` subprogram
...........................
The first step in this |metagene| analysis is to examine a genome annotation,
find all of the features of interest (transcripts), and align them at some
interesting sub-feature (a start codon). This is handled by the |metagene|
``generate`` subprogram. Specifically, |metagene| ``generate``:

 #. Groups transcripts by gene, so each start codon is counted once

 #. If a gene has multiple start codons, that gene is excluded.
    
    If a gene has a single start codon, then |metagene| finds the *maximal
    spanning window* of the gene, defined as the largest possible
    window surrounding the start codon in which all transcripts from
    that gene map to identical genomic positions. This prevents any
    positional ambiguity from entering the average.

    .. TODO : make image for this

 #. Once maximal spanning windows are defined for each gene, |metagene|
    determines the location of the start codon relative to the new window
 
 #. Saves this info to a `BED`_ file for genome browsing, and a text file
    (called an *ROI file*) for further processing.

The ``generate`` program only needs to be run once per landmark of interest
(i.e. once for stop codons, once for start codons), and only needs to be 
re-run if the genome :term:`annotation` changes (e.g. due to revisions,
additions, or deletions of gene/transcript models).

The ``generate`` program may be called from the terminal:

 .. code-block:: shell

    $ metagene generate yeast_cds_start --landmark cds_start --annotation_files sgd_plus_utrs_chrI.gtf


The ``count`` subprogram
........................
Once ``generate`` has made an ROI file, |metagene|'s ``count`` subprogram can be 
used to make metagene averages. Specifically, ``count`` takes the following
steps:

 #. For each ROI in the ROI file:

     #. fetch a vector of counts at each position from the sample dataset

     #. If the vector of counts has enough counts within a user-specified
        *normalization region*, include it. Otherwise, exclude the vector
     
     #. Normalize the vector by the number of counts in the *normalization region*

 #. Construct a metagene average by taking the median over all
    normalized vectors at each position, excluding at any specific position
    all vectors that do not cover that position (because the maximal spanning
    window for that gene was too small).

 #. Save the metagene average and the number of genes counted at each position
    to a tab-delimited text file.

To call the ``count`` program, type into a terminal window. In this example
``--threeprime --offset 15`` specify our :term:`mapping rule` for the
:term:`P-site offset`. 

 .. code-block:: shell

    $ metagene count yeast_cds_start_rois.txt SRR1562907 --count_files SRR1562907_chrI.bam --threeprime --offset 15


A number of files are created. ``SRR1562907_metagene_profile.txt`` is the
final product. It contains three columns:

  #. *x:* an X-coordinate indicating the distance in nucleotides from
     the start codon
  
  #. *metagene_profile:* the value of the metagene average at *x*
       
  #. *regions_counted:* the number of regions included in the average at *x*


The ``chart`` subprogram
........................
For convenience, a chart subprogram is included. It can plot multiple
metagene profiles (each from a run of the ``count`` subprogram) in
a single plot:

 .. code-block:: shell
 
    $ metagene chart SRR1562907_cds_start.png SRR1562907_metagene_profile.txt --landmark "Start codon" --title "Metagene chart title"

This produces the image:

    [TODO]: insert image



.. _concepts-metagene-roll-your-own:

Beyond start and stop codons: defining your own window functions
----------------------------------------------------------------

:term:`Metagene averages <metagene average>` can be useful for other questions,
regions, and experimental data types. For this reason, |metagene| offers tools
to create maximal spanning windows surrounding any feature of interest.

To do so,
it requires a *window function* that can identify and build a window around
the feature of interest in each individual region (for example, a transcript).
If regions define a `"gene_id"` attribute in their `attr` dictionaries,
then their windows grouped by shared values for `"gene_id"` when their maximal
spanning windows are made.

|metagene| comes with two window functions, which its command-line program uses:

  - :func:`~yeti.bin.metagene.window_cds_start`, for defining windows
    surrounding start codons

  - :func:`~yeti.bin.metagene.window_cds_stop`, for defining windows
    surrounding stop codons

However, it is possible to design window functions that focus on any arbitrary 
feature of interest, such as:

  - a specific nucleic acid sequence

  - features (e.g. spikes, valleys, step changes in density) in
    :term:`high-throughput sequencing` data or other quantitative 
    data

  - SNPs or other mutations

To perform a metagene analysis around such a feature, you need
to write a window function that identifies that *feature* within your
*regions of interest* (in the start codon example, a *start codon* within
each *transcript*).

Then, the function :func:`yeti.bin.metagene.do_generate` can use your
window function to generate maximal spanning windows.

Window functions must take the following parameters, in order:

    `roi` : |SegmentChain|
        Input ROI for which to generate a window.
        If `"gene_id"` is defined in `roi.attr`, then the 
        all windows from all ROIs sharing the same `"gene_id"`
        as this one will be used to generate a single
        maximal spanning window covering all of them.

    `flank_upstream` : ``int``
        Nucleotide length upstream of the feature of interest
        to include in the maximal spanning window , if `roi` has
        such a feature

    `flank_downstream` : ``int``
        Nucleotide length downstream of the feature of interest
        to include in the maximal spanning window, if `roi` has
        such a feature
    
    `ref_delta` : ``int``, optional
        Offset in nucleotides from the feature of interest to 
        the reference point at which all maximal spanning window
        count vectors will be aligned when the metagene average
        is calculated. If `0`, the feature of interest is the
        reference point. (Default: `0`)


Window functions must return the following values, in order:

    |SegmentChain|
        Window surrounding feature of interest if `roi` has such a feature.
        Otherwise, return a zero-length |SegmentChain| 
    
    ``int``
        offset to align window with all other windows, if `roi` itself
        wasn't long enough in the 5\' direction to include the entire
        distance specified by `flank_upstream`. Otherwise `0`. 

    (``str``, ``int``, ``str``)
        Genomic coordinate of reference point as `(chromosome name, coordinate, strand)`


Here is a window function that produces windows surrounding transcription
start sites::

    >>> def window_transcript_start(roi,flank_upstream,flank_downstream,ref_delta=0):
    >>> """Window function for metagenes surrounding roiion start sites
    >>>
    >>> Returns
    >>> -------
    >>> SegmentChain
    >>>     Window surrounding transcript start site
    >>>
    >>> int
    >>>     Offset to align window with all other windows
    >>>
    >>> (str,int,str)
    >>>     Genomic coordinate of transcription start site as *(chromosome name, coordinate, strand)*
    >>> """
    >>> ref_point = roi.get_genomic_coordinate(0)
    >>> segs = roi.get_subchain(0,flank_downstream)
    >>>
    >>> if ref_point[2] == "+":
    >>>     new_segment_start = ref_point[1] - flank_upstream
    >>>     new_segment_end = roi.get_genomic_coordinate(flank_downstream)
    >>>     offset = 0
    >>> else:
    >>>     new_segment_start = roi.get_genomic_coordinate(flank_downstream)
    >>>     new_segment_end = ref_point[1] + flank_upstream
    >>>     if roi.get_length() < flank_downstream:
    >>>         offset = flank_downstream - roi.get_length()
    >>>
    >>> outside_segment = GenomicSegment(ref_point[0],
    >>>                                  new_segment_start,
    >>>                                  new_segment_end,
    >>>                                  ref_point[2])
    >>> segs.append(outside_segment)
    >>> new_chain = SegmentChain(*tuple(segs))
    >>>
    >>> return new_chain, offset, ref_point

.. TODO : test this window function

Here is a window function that produces windows surrounding the highest spike
in read density in a transcript. Note, it uses data structures in the global
scope::

    >>> TODO


Once your window function is written, you can generate maximal spanning windows::


    >>> from yeti.bin.metagene import do_generate

    >>> # include 100 nucleotides up- and downstream of feature
    >>> flank_upstream = flank_downstream = 100

    >>> data_table, segment_chains = do_generate(transcripts,mask_hash,
    >>>                                          flank_upstream,flank_downstream,
    >>>                                          landmark_func=window_largest_spike)

:meth:`~yeti.bin.metagene.do_generate` returns an |ArrayTable| (similar to a :class:`pandas.DataFrame`)
of data and a list of |SegmentChains| corresponding to the maximal spanning windows for each row
in the |ArrayTable|. It is useful to save these in the same formats that the |metagene| ``generate``
program uses::

    >>> with open("my_roi_file.txt","w") as roi_fh:
    >>>     data_table.to_file(roi_fh)
    >>>     data_fh.close()
    >>>     
    >>> with open("my_rois.bed","w") as bed_fh:
    >>>     for roi in segment_chains:
    >>>         bed_fh.write(roi.as_bed())
    >>> 
    >>>     bed_fh.close()
 
Then, you can use the |metagene| ``count`` subprogram to do your counting:

 .. code-block:: shell

     $ metagene count my_roi_file.txt SRR1562907_custom_metagene --count_files SRR1562907_chrI.bam --threeprime --offset 15



See also
--------
  - Module documentation for |metagene| program

    
