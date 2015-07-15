Tour
====

This document contains a brief overview of the :ref:`tour-scripts` and
:ref:`object types <tour-data-structures>` included in :data:`yeti`. More detail
on all of these may be found in the :doc:`module documentation <generated/yeti>`.


.. _tour-scripts:

Command-line scripts
--------------------

:data:`yeti` includes a handful of :mod:`scripts <yeti.bin>` that perform common
:term:`high-throughput sequencing` and :term:`ribosome profiling` analyses.
These include:

  - :doc:`Measuring read density </generated/yeti.bin.cs>` in regions
    of interest (for example, to obtain gene expression values)

  - :doc:`Determining P-site offsets </generated/yeti.bin.psite>` for
    :term:`ribosome profiling` experiments

  - :doc:`Performing metagene analyses </generated/yeti.bin.metagene>` for
    any sort of :term:`high-throughput sequencing` experiment

  - :doc:`Creating browser tracks </generated/yeti.bin.make_wiggle>` 
    from `BAM`_ or `bowtie`_ files, after applying optional read :term:`mapping rules`
    (e.g. for P-site assignment, in the case of ribosome profiling data) 

  - and others!

For a complete list, as well as examples, see the :mod:`command-line script documentation <yeti.bin>`



.. _tour-data-structures:

Classes & data structures
-------------------------

:data:`yeti` defines a number of classes that facilitate sequencing
analyses. These are:

    ======================================================     =========================================
    **Class**                                                  **Purpose**
    ------------------------------------------------------     -----------------------------------------
    :ref:`tour-genomic-segment`, :ref:`tour-segment-chain`     Represent genomic features (e.g. mRNAs, genes, SNPs, stop codons) as Python objects

    :ref:`GenomeArray <tour-genome-array>` & subclasses        Array-like object that maps quantitative values (e.g. read counts, phylogenetic conservation)
                                                               to corresponding genomic coordinates.

    :ref:`GenomeHash <tour-genome-hash>` & subclasses          Array-like object that indexes genomic features by genomic coordinates, 
                                                               for quick lookup of features that overlap or cover a region.
    ======================================================     =========================================

In the examples below, we'll be using a small :doc:`test_dataset` for yeast chromosome I.

-------------------------------------------------------------------------------

.. _tour-genomic-segment:

|GenomicSegment|
................
|GenomicSegments| are the fundamental building block of genomic features.
They are defined by a chromosome name, a start coordinate, and end coordinate,
and a strand. On their own, they are not that interesting. However, they
can be used to build :ref:`segment-chain`, which are interesting.


.. _tour-segment-chain:

|SegmentChains| & |Transcripts|
...............................

|SegmentChain| & its subclass |Transcript| model genomic features. They are
constructed from zero or more |GenomicSegments|, and therefore can represent
even discontinuous genomic features, such as transcripts or gapped alignments,
in addition to continuous features (e.g. single exons).
	
|SegmentChain| and its subclasses provide methods for:
	
  - converting coordinates between the genome and the spliced space of the
    |SegmentChain|

  - fetching genomic sequence, read alignments, or count data over
    the |SegmentChain|, in its own 5' to 3' direction, automatically
    accounting for splicing of the segments and, for reverse-strand
    features, reverse-complementing the sequence

  - slicing or fetching sub-regions of a |SegmentChain|
      
  - testing for equality, inequality, overlap, containment, or coverage
    of other |SegmentChain| or |GenomicSegment| objects

  - exporting to `BED`_, `GTF2`_, or `GFF3`_ formats, for use with other
    software packages or within a genome browser

|SegmentChains| and |Transcripts| can be constructed manually from zero or more
|GenomicSegments| and any optional keywords, which will be stored in the
|SegmentChain|'s `attr` dictionary::

    >>> from yeti.genomics.roitools import *
    >>> exon1 = GenomicSegment("chrI",129237,130487,"+")
    >>> exon2 = GenomicSegment("chrI",130531,130572,"+")
    >>> SegmentChain(exon1,exon2,ID="YAL013W",alias="DEP1")
    <SegmentChain intervals=2 bounds=chrI:129237-130572(+) name=YAL013W>

    >>> dep1 = Transcript(exon1,exon2,ID="YAL013W",alias="DEP1",cds_genome_start=129270,cds_genome_end=130484)
    >>> dep1
    <Transcript intervals=2 bounds=chrI:129237-130572(+) name=YAL013W>
    
    >>> dep1.attr
    {'ID': 'YAL013W',
     'alias': 'DEP1',
     'cds_genome_end': 130484,
     'cds_genome_start': 129270,
     'type': 'mRNA'}


More often, |SegmentChains| and |Transcripts| are loaded from :term:`annotation` files
(see :mod:`yeti.readers`, which parsers the :term:`annotation` files). To assemble
transcripts from exons and coding regions in a `GTF2`_ file::
 
    >>> from yeti.readers.gff import GTF2_TranscriptAssembler

    >>> # get an iterator over transcripts in file
    >>> reader = GTF2_TranscriptAssembler(open("sgd_plus_utrs_chrI.gtf"))

    >>> # do something with transcripts. here we just look at their names & attribute dictionaries
    >>> for transcript in reader:
            print(transcript.get_name() + ":\t" + str(transcript.attr))
    YAL069W_mRNA:	{'cds_genome_end': 646, 'name': 'YAL069W', 'gene_id': 'YAL069W', 'utr5_source': 'estimated', 'source': '.', 'transcript_id': 'YAL069W_mRNA', 'cds_genome_start': 334, 'phase': '.', 'utr3_source': 'estimated', 'gene_aliases': 'YAL069W', 'score': '.', 'type': 'mRNA', 'ID': 'YAL069W_mRNA'}
    YAL068W-A_mRNA:	{'cds_genome_end': 789, 'name': 'YAL068W-A', 'gene_id': 'YAL068W-A', 'utr5_source': 'estimated', 'source': '.', 'transcript_id': 'YAL068W-A_mRNA', 'cds_genome_start': 537, 'phase': '.', 'utr3_source': 'estimated', 'gene_aliases': 'YAL068W-A', 'score': '.', 'type': 'mRNA', 'ID': 'YAL068W-A_mRNA'}
    YAL068C_mRNA:	{'cds_genome_end': 2169, 'name': 'PAU8', 'gene_id': 'YAL068C', 'utr5_source': 'estimated', 'source': '.', 'transcript_id': 'YAL068C_mRNA', 'cds_genome_start': 1809, 'phase': '.', 'utr3_source': 'estimated', 'gene_aliases': 'PAU8,seripauperin PAU8', 'score': '.', 'type': 'mRNA', 'ID': 'YAL068C_mRNA'}


|SegmentChains| and |Transcripts| can convert coordinates between the transcript
and the genome::

    >>> # load transcripts into a dictionary keyed on transcript ID
    >>> transcript_dict = { X.get_name() : X for X in GTF2_TranscriptAssembler(open("sgd_plus_utrs_chrI.gtf")) }

    >>> # we'll use the two-exon, minus-strand gene TFC3 as an example
    >>> tfc3 = transcript_dict["YAL001C_mRNA"]
    >>> tfc3
    <Transcript intervals=2 bounds=chrI:147529-151186(-) name=YAL001C_mRNA>

    >>> # get genomic coordinate of 89th nucleotide from 5' end of TFC3
    >>> # right before the splice junction
    >>> tfc3.get_genomic_coordinate(89)
    ('chrI', 151096, '-')
    
    >>> # get genomic coordinate of 90th nucleotide from 5' end of TFC3
    >>> # right after the splice junction
    >>> tfc3.get_genomic_coordinate(90)
    ('chrI', 151005, '-')

    >>> # and the inverse operation also works
    >>> tfc3.get_segmentchain_coordinate('chrI', 151005, '-')
    90

.. _tour-get-counts:

One of the most convenient things |SegmentChains| can do is to fetch vectors of
data covering each position in the chain from the 5' to 3' end (data is held
in |GenomeArrays| which are explained :ref:`below <tour-genome-array>`). 
In this example, we count how many 5' ends of sequencing reads appear at each
position in the chain::

    >>> from yeti.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory
    >>> import pysam

    >>> # load read alignments, and map them to 5' ends
    >>> alignments = BAMGenomeArray([pysam.Samfile("SRR1562907_chrI.bam","rb")])
    >>> alignments.set_mapping(FivePrimeMapFactory())

    >>> # fetch the number of 5' ends of alignments at positions 300-320
    >>> tfc3.get_counts(alignments)[300:320]
    array([ 0.,  0.,  0.,  2.,  0.,  0.,  1.,  0.,  1.,  0.,  0.,  1.,  1.,
            0.,  0.,  0.,  0.,  0.,  0.,  0.])


It is also possible to fetch sub-sections of a |Transcript| or |SegmentChain|
as a new |SegmentChain|::

    >>> # take first 200 nucleotides of TFC3 mRNA
    >>> subchain = tfc3.get_subchain(0,200)
    >>> subchain
    <SegmentChain intervals=2 bounds=chrI:150896-151186(-) name=YAL001C_mRNA>

|Transcript| includes several convenience methods to fetch 5' UTRs, coding regions,
and 3'UTRs from coding transcripts::

    >>> tfc3.get_utr5()
    <SegmentChain intervals=1 bounds=chrI:151166-151186(-) name=YAL001C_mRNA>

    >>> tfc3.get_cds()
    <SegmentChain intervals=2 bounds=chrI:147596-151166(-) name=YAL001C_mRNA>


|SegmentChain| and its subclasses can also fetch their sequences from dictionaries
of strings or :class:`bio.SeqRecord.SeqRecord` objects. These sequences will
automatically be spliced and reverse-complemented, as necessary::

    >>> from Bio import SeqIO
    >>> genome = SeqIO.to_dict(SeqIO.parse(open("chrI.fa"),"fasta"))
    >>> tfc3.get_cds().get_sequence(genome)
    'ATGGTACTGACGATTTATCCTGACGAACTCGTACAAATAGTGTCTGATAAAATTGCTTCAAATAAGGGAAAAATCACTTTGAATCAGCTGTGGGATATATCTGGTAAATATT
    # rest of output omitted


|SegmentChains| and |Transcripts| can do a lot more. For a complete description
of their functions, attributes, et c, see the documentation for |SegmentChain|
and |Transcript| in the :py:mod:`~yeti.genomics.roitools` module.
    
-------------------------------------------------------------------------------

.. _tour-genome-array:

|GenomeArray| & its subclasses
..............................
|GenomeArrays| store count data at each position in the genome. Data can be
imported from count files (e.g. `Wiggle`_, `bedGraph`_) as well as alignment files
(in `bowtie`_ or `BAM`_ format). For very large genomes a sparse implementation
is provided by |SparseGenomeArray|. A |BAMGenomeArray| is provided for
:term:`read alignments` in `BAM`_ format.

When importing :term:`read alignment`, users can specify a :term:`mapping functions <mapping function>`
to determine how the alignments should be converted into counts (e.g., to their
fiveprime ends, threeprime ends, or, something more complex). :data:`yeti` already
includes a number of mapping functions to map read alignments:

  - to their fiveprime or threeprime ends, with or without offsets from
    the end of the read. These offests can be constant, or a function of 
    read length (e.g. for :term:`P-site mapping` for :term:`ribosome profiling data`). 
     
  - fractionally over their entire lengths (e.g. for RNA-seq)
   
  - fractionally to all positions covered by a central portion of the read
    alignment, after excluding a user-defined number of positions on each
    send of the read (as in ribosome profiling data from *E. coli*
    :cite:`Oh2011` or *D. melanogaster* :cite:`Dunn2013`).


As we saw :ref:`above <tour-get-counts>`, |GenomeArrays| are most often
called  using the :meth:`~yeti.genomics.roitools.SegmentChain.get_counts`
method of |SegmentChains| or |Transcripts|. But, they have a few other abilities,
too. For example, :term:`mapping functions <mapping function>` for 
|BAMGenomeArrays| can be changed at runtime::

    >>> from yeti.genomics.genome_array import FivePrimeMapFactory, ThreePrimeMapFactory
    
    >>> alignments.set_mapping(FivePrimeMapFactory())
    >>> tfc3.get_cds_().get_counts(alignments)[:50]
    array([ 3.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,
            0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,
            1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    >>> # change to mapping with 15 nucleotide offset from 5' end
    >>> alignments.set_mapping(FivePrimeMapFactory(offset=15))
    >>> tfc3.get_cds_().get_counts(alignments)[:50]
    array([ 0.,  0.,  3.,  2.,  1.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  3.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

    >>> # change to mapping from 3' end, with no offset
    >>> alignments.set_mapping(ThreePrimeMapFactory())
    >>> tfc3.get_cds().get_counts(alignments)[:50]
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  8.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  2.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])


|GenomeArrays| can be exported to `wiggle`_ or `bedGraph`_ files::

    >>> # export minus strand as a bedgraph file
    >>> with open("alignments_rc.wig","w") as fout:
            alignments.to_bedgraph(fout,"my_trackname","-")


And `wiggle`_ or `bedGraph`_ files can be reimported. Both use the
:meth:`~yeti.genomics.genome_array.GenomeArray.add_from_wiggle` method::

    >>> # reimport plus strand data
    >>> new_data = GenomeArray()
    >>> new_data.add_from_wiggle(open("alignments_rc.wig"),"-")
    
    >>> # should be the same as using 3' mapping, because this was
    >>> # the mapping rule used during export
    >>> tfc3.get_cds().get_counts(new_data)[:50]
    array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  8.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  2.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])


For further information, see:

  - The module documentation for :py:mod:`~yeti.genomics.genome_array`

  - In-depth discussion of :doc:`mapping rules <concepts/mapping_rules>`

-------------------------------------------------------------------------------

.. _tour-genome-hash:

|GenomeHash|, |BigBedGenomeHash|, and |TabixGenomeHash|
.......................................................

Often one needs to know whether any features overlap a specific region in the
genome, for example, to find transcripts that overlap one another.

However, it would be inefficient to scan an entire file to find the
overlapping features, or to test whether two features overlap if we already
know they are too far apart in the genome.

|GenomeHash| and its subclasses avoid this problem by indexing features
by location. A |GenomeHash| may be created from a list or dictionary of features
(e.g. |SegmentChains| or |Transcripts|) in memory, or directly loaded from a
genome annotation (in `BED`_, `GTF2`_, `GFF3`_, or `PSL`_ format)::

    >>> from yeti.genomics.genome_hash import GenomeHash 
    >>> my_hash = GenomeHash(transcript_dict)

Having made a |GenomeHash|, we can ask what is where in the genome. For
example, to find all features between bases 10000-20000 on the plus
strand of chromosome *chrI*::

    >>> from yeti.genomics.genome_hash import GenomeHash 
    >>> my_hash = GenomeHash(transcript_dict)
    
    >>> roi = GenomicSegment("chrI",10000,20000,"+")
    >>> my_hash[roi]
    [<Transcript intervals=1 bounds=chrI:9979-10540(+) name=YAL066W_mRNA>,
     <Transcript intervals=1 bounds=chrI:11934-12567(+) name=YAL064W-B_mRNA>]

Or on both strands::

    >>> my_hash.get_overlapping_features(roi,stranded=False)
    [<Transcript intervals=1 bounds=chrI:11423-12062(-) name=YAL065C_mRNA>,
     <Transcript intervals=1 bounds=chrI:9979-10540(+) name=YAL066W_mRNA>,
     <Transcript intervals=1 bounds=chrI:11934-12567(+) name=YAL064W-B_mRNA>,
     <Transcript intervals=1 bounds=chrI:13221-13854(-) name=YAL064C-A_mRNA>]
    
Does anything interesting overlap *TFC3*?

 .. code-block:: python

    >>> my_hash[tfc3]
    [<Transcript intervals=2 bounds=chrI:147529-151186(-) name=YAL001C_mRNA>]
    # nope, just TFC3 itself.

For more information, see the module documentation for :mod:`~yeti.genomics.genome_hash`.

-------------------------------------------------------------------------------


See also
--------
For an in-depth discussion of these data structures, see:
	
  - :doc:`Example <examples>` analyses

  - Detailed :ref:`module documentation <modindex>`, especially:

      - :mod:`yeti.genomics.roitools`

      - :mod:`yeti.genomics.genome_array`

      - :mod:`yeti.genomics.genome_hash`

      - :mod:`yeti.readers`
