Read phasing in :term:`ribosome profiling`
==========================================

During each cycle of peptide elongation, the ribosome takes a 3-nucleotide
step along its mRNA. This physical process creates
:term:`triplet periodicity` in :term:`ribosome profiling` data, which
becomes visible when the :term:`read alignments`  are mapped
to their :term:`P-site offsets <P-site offset>` (:cite:`Ingolia2009`):

 .. TODO: phasing figure

 .. figure: 
    :alt: Ribosome phasing genome browser examples
    :figclass: captionfigure

    :term:`Triplet periodicity` across a coding region in :doc:`/test_dataset`


:term:`Triplet periodicity` or :term:`sub-codon phasing`, is a useful
feature, becuase it can be used to infer the reading frame(s) in which
a coding region is decoded:

 .. TODO: insert phasing chart figure

 .. figure:
    :alt: Phasing differs between reading frames
    :figclass: captionfgure

    :term:`triplet periodicity` provides unique signatures of reading frames


 .. note::

    Heavily biased nucleases -- for example, micrococcal
    nuclease, used for :term:`ribosome profiling` of  *E. coli* (:cite:`Oh2011`)
    and *D. melanogaster* (:cite:`Dunn2013`) -- introduce uncertainty into
    P-site mapping, which obscures :term:`phasing <sub-codon phasing>`
    in these datasets.


In this tutorial, we show how to measure :term:`sub-codon phasing`

  - :ref:`Manually <examples-phasing-manual>`, in an interactive Python session
  - :ref:`From the command-line <exampels-phasing-script>`, using the
    |phase_by_size| script


The examples below use the :doc:`/test_dataset`.


 .. _examples-phasing-manual

Examining phasing manually
..........................

In this section, we tabulate :term:`sub-codon phasing` manually.
First, we need to open the read alignments and transcript annotation:

 .. code-block:: python

    >>> import pysam
    >>> from plastid.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory
    >>> from plastid.readers.bed import BED_Reader

    >>> # retrieve an iterator over transcripts
    >>> transcripts = BED_Reader(open("merlin_orfs.bed"),return_type=Transcript)

    >>> # open read alignments and map to P-sites
    >>> alignments = BAMGenomeArray([pysam.Samfile("SRR609197_riboprofile.bam","rb")])
    >>> alignments.set_mapping(FivePrimeMapFactory(offset=14))

:term:`Ribosome-protected footprints <footprint>` of varying lengths exhibit variable
phasing. So, we'll look at the most highly-phased population of reads, 33-mers. 
To do so, we'll add a size filter:

    >>> from plastid.genomics.genome_array import SizeFilterFactory
    >>> size_filter = SizeFilterFactory(min=33,max=34)
    >>> >>> alignments.add_filter("size",size_filter)


Next, we can count phasing:

 .. code-block:: python

    >>> # create a holder for phasing
    >>> phasing = numpy.zeros(3)
    
    >>> # start codons are hyper-phased; stop codons can have differnt
    >>> # phasing or even be de-phased depending on experimental protocol
    >>> # so, we'll ignore 5 codons after the start, and 5 before the stop
    >>> codon_buffer = 5*3

    >>> # count
    >>> for my_transcript in transcripts:
    >>>     cds = my_transcript.get_cds()
    >>>     # if transcript is coding
    >>>     if len(cds) > 0: 
    >>>         try:
    >>>
    >>>             # get numpy.ndarray of counts in coding region
    >>>             counts = cds.get_counts(alignments)[codon_buffer:-codon_buffer]
    >>>
    >>>             # reshape to Nx3, where N = number of codons
    >>>             counts = counts.reshape((len(counts)/3,3))
    >>>
    >>>             # sum over codon positions to get a 3-vector,
    >>>             # and add to data holder
    >>>             phasing += counts.sum(0)
    >>>
    >>>         except: # raise exception if coding region is not n*3 nucleotides long
    >>>             print("Length (%s nt) of CDS for `%s` contains partial codons. Frameshift?" % (len(counts),my_transcript.get_name()))

    >>> # compute fraction of phased reads
    >>> phasing_proportions = phasing.astype(float) / phasing.sum()
    >>> phasing_proportions
    array([ 0.51042163,  0.29362327,  0.19595509])

 .. note::

    If the transcript annotation includes multiple transcript isoforms
    for the same gene, codons that appear in more than one isoform will
    be double-counted in the phasing estimate. This may be avoided by
    filtering the annotation file ahead of time.
    
    If the annotation file contains overlapping coding regions which appear
    in different frames, including these in the phasing tabulation will 
    under-estimate phasing. It makes sense to exclude such areas using a
    :term:`mask file`.


 .. _examples-phasing-script

Measuring :term:`phasing <sub-codon phasing>` using the |phase_by_size| script
..............................................................................

The |phase_by_size| script automates the calculations described in 
:ref:`examples-phasing-manual`, calculating phasing separately for
:term:`read alignments` of each length.

The command line below examines phasing in 
the :term:`ribosome profiling` dataset ``SRR609197_riboprofile.bam``,
estimating the P-site as 14 nucleotides from the 5' end of each read.
In addition, we exclude 5 codons near the start and stop codons (
via the ``--codon_buffer`` argument), because these are often hyper-phased
compared to coding regions:

  .. code-block:: shell

     $ phase_by_size SRR609197 \
                     --count_files SRR609197_riboprofile.bam \
                     --annotation_files merlin_orfs.bed \
                     --annotation_format BED \
                     --fiveprime --offset 14 \
                     --codon_buffer 5 \
                     --min_length 25 --max_length 35

|phase_by_size| will create a text file showing the proportion of reads
whose P-sites map to each codon position for each read length (columns
*phase0, phase1,* & *phase0*) and the proportion of total reads that
each read length represents (column *fraction_reads_counted*):

 .. code-block:: shell

    #read_length    reads_counted    fraction_reads_counted    phase0      phase1      phase2
    25              6511             0.009640                  0.326832    0.327599    0.345569
    26              9952             0.014735                  0.385953    0.295217    0.318830
    27              17636            0.026111                  0.320934    0.282717    0.396348
    28              42976            0.063629                  0.251792    0.381794    0.366414
    29              93754            0.138809                  0.309309    0.370971    0.319720
    30              148400           0.219716                  0.318733    0.367635    0.313632
    31              155684           0.230501                  0.336624    0.421713    0.241663
    32              118565           0.175543                  0.445578    0.374141    0.180281
    33              58761            0.087000                  0.511121    0.299076    0.189803
    34              18818            0.027861                  0.508237    0.276597    0.215166
    35              4360             0.006455                  0.514450    0.236468    0.249083


 .. note::

    At present, |phase_by_size| only supports :term:`read alignments`
    in `BAM`_ format.

-------------------------------------------------------------------------------

See also
--------
  - :doc:`/examples/p_site` for a discussion on how to determine the
    :term:`P-site offsets <P-site offset>` to use for a given
    :term:`ribosome profiling` dataset.
  - :doc:`/concepts/mapping_rules` for a discussion on how to apply
    :term:`P-site offsets <P-site offset>` or other mapping rules
