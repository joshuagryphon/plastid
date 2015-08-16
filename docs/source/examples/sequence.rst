Manipulating genomic sequence
=============================
In this tutorial we will fetch the sequences of genomic `features <feature>`_,
using the :doc:`/test_dataset`.


Fetching sequence from |SegmentChains|
--------------------------------------
|SegmentChains| and |Transcripts| can fetch their own genomic sequence via
their :meth:`~yeti.genomics.roitools.SegmentChain.get_sequence` method.
Sequences are fully-spliced, if the |SegmentChain| contains one or more
|GenomicSegment|::

    >>>

For reverse-strand features, sequences are automatically reverse-complemented::

    >>>

For convenience, |SegmentChain| and |Transcript| also implement
:meth:`~yeti.genomics.roitools.SegmentChain.get_fasta`, which formats output
in `FASTA`_ format::

    >>>


Supported sequence file formats
-------------------------------
Genomic sequence can be stored in a number of formats. Two of the most common are:

  - `FASTA`_: a non-indexed, text-base file of one or more sequences.
  - `2bit <twobit>`_: an indexed, binary format created for large (e.g. mammalian)
     genomes. Beacause `2bit <twobit>`_ files are binary and indexed, they use
     far less memory than `FASTA`_ files.

:data:`yeti` doesn't implement any of its own file readers, but it is compatible
with any parser that returns sequences as a dictionary-like object of string-like
objects. This allows users to use implementations in `BioPython`_ (for `FASTA`_, 
EMBL, & GenBank formats) and `twobitreader`_ (for `2bit <twobit>`_) files.

The objects returned from these readers can also be passed to
:meth:`~yeti.genomics.roitools.SegmentChain.get_sequence`::

    >>> # twobit reader example
    >>>
    >>>

    >>> # BioPython example
    >>>


Manipulating sequence
---------------------
Tools for further manipulating sequence (e.g. reverse-complementing, translating)
are supplied in `BioPython`_'s `Seq`_ and `SeqRecord`_ objects::

    >>> # SeqRecord examples
    >>>


-------------------------------------------------------------------------------

See also
--------
  - `BioPython`_ documentation for manipulation of nucleic acid sequences
  - `twobitreader`_ documentation
  - `UCSC file format FAQ`_ for details on sequence file formats
