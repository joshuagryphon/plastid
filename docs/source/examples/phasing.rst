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
    and *D. melanogaster* (:cite:`Dunn2014`) -- introduce uncertainty into
    P-site mapping, which obscures :term:`phasing <sub-codon phasing>`
    in these datasets.


In this tutorial, we show how to measure :term:`sub-codon phasing`

  - :ref:o`Manually <examples-phasing-manual>`, in an interactive Python session
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
    >>> from yeti.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory
    >>> from yeti.readers.bed import BED_Reader

    >>> # open read alignments and map to P-sites
    >>> alignments = BAMGenomeArray([pysam.Samfile("","rb")])
    >>> alignments.set_mapping(FivePrimeMapFactory(offset=14))

    >>> retrieve an iterator over transcripts
    >>> transcripts = BED_Reader(open("merlin_orfs.bed"),return_type=Transcript)

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
:ref:`examples-phasing-manual`, calculating phasing seprately for
:term:`read alignments` of each length.

The command line below examines phasing in 
the :term:`ribosome profiling` dataset ``SRR609197_riboprofile.bam``,
estimating the P-site as 14 nucleotides from the 5' end of each read.
In addition, we exclude 5 codons near the start and stop codons (
via the ``--codon_buffer`` argument), because these are often hyper-phased
compared to coding regions:

  .. code-block:: shell

     $ phase_by_size SRR609197_phasing \
                     --count_files SRR609197_riboprofile.bam \
                     --annotation_files merlin_orfs.bed \
                     --annotation_file_format BED \
                     --fiveprime --offset 14 \
                     --codon_buffer 5

|phase_by_size| will create a text file showing the raw number
and proportion of reads whose P-sites map to each codon position
for each read length, an image that visually shows the same data:

 .. figure:: TODO
    :figclass:captionfigure
    :alt: Output of |phase_by_size| script

    Sample of graphical output of |phase_by_size|

 .. TODO: remove note if otehr file format support added to phase_by_size

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
