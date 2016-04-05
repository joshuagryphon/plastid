Tour
====

This document contains a brief overview of the :ref:`tour-scripts` and
:ref:`object types <tour-data-structures>` included in :data:`plastid`. Complete
documentation may be found in :doc:`module documentation <generated/plastid>`.

.. contents::
   :local:



.. _tour-scripts:

Command-line scripts
--------------------

:data:`plastid` includes a handful of :mod:`scripts <plastid.bin>` for common
:term:`high-throughput sequencing` and :term:`ribosome profiling` analyses.
These include, among others:

 - :doc:`Measuring read density </generated/plastid.bin.cs>` in regions
   of interest (for example, to obtain gene expression values)

 - :doc:`Determining P-site offsets </generated/plastid.bin.psite>` for
   :term:`ribosome profiling` experiments

 - Performing :doc:`metagene analyses </generated/plastid.bin.metagene>`

 - :doc:`Creating browser tracks </generated/plastid.bin.make_wiggle>` 
   from `BAM`_ or `bowtie`_ files, after applying read :term:`mapping rules`
   to transform the alignments (e.g. for P-site assignment of 
   :term:`ribosome profiling` data) 

For a complete list, see :doc:`/scriptlist`.



.. _tour-data-structures:

Classes & data structures
-------------------------

:data:`plastid` defines a number of classes to facilitate sequencing analyses:

   =======================================================    ===============================================
   **Class**                                                  **Purpose**
   -------------------------------------------------------    -----------------------------------------------
   :ref:`tour-genomic-segment`, :ref:`tour-segment-chain`     Represent genomic :term:`features <feature>`
                                                              (e.g. mRNAs, genes, SNPs, stop codons) as
                                                              Python objects

   :ref:`GenomeArray <tour-genome-array>` & its subclasses    Map quantitative values or
                                                              :term:`read alignments` to genomic coordinates.

   :ref:`GenomeHash <tour-genome-hash>` & its  subclasses     Index genomic :term:`features <feature>` by
                                                              genomic coordinates, for quick lookup of
                                                              :term:`features <feature>` that overlap or
                                                              cover a region.
   =======================================================    ===============================================

In the examples below, we'll be using a small :doc:`test_dataset` covering the human cytomegalovirus (hCMV) genome (:cite:`Stern-Ginossar2012`).

-------------------------------------------------------------------------------

.. _tour-genomic-segment:

|GenomicSegment|
................
|GenomicSegments| are the fundamental building block of genomic
:term:`features <feature>`. They are defined by:

 - a chromosome name
  
 - a start coordinate
  
 - an end coordinate
  
 - a strand:
  
    - '+' for forward-strand features
    - '-' for reverse-strand features
    - '.' for unstranded features
      
On their own, |GenomicSegments| are not very interesting. However, they
can be used to build :ref:`SegmentChains <tour-segment-chain>`, which are interesting.

-------------------------------------------------------------------------------

.. _tour-segment-chain:

|SegmentChain| & |Transcript|
.............................

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

   >>> from plastid import GenomicSegment, SegmentChain, Transcript
   >>> exon1 = GenomicSegment("chrI",129237,130487,"+")
   >>> exon2 = GenomicSegment("chrI",130531,130572,"+")
   >>> SegmentChain(exon1,exon2,ID="YAL013W",alias="DEP1")
   <SegmentChain segments=2 bounds=chrI:129237-130572(+) name=YAL013W>

   >>> dep1 = Transcript(exon1,exon2,ID="YAL013W",alias="DEP1",cds_genome_start=129270,cds_genome_end=130484)
   >>> dep1
   <Transcript segments=2 bounds=chrI:129237-130572(+) name=YAL013W>
    
   >>> dep1.attr
   {'ID': 'YAL013W',
    'alias': 'DEP1',
    'cds_genome_end': 130484,
    'cds_genome_start': 129270,
    'type': 'mRNA'}


More often, |SegmentChains| and |Transcripts| are loaded from :term:`annotation`
files (see :mod:`plastid.readers`)::
 
   >>> from plastid import BED_Reader

   >>> # get an iterator over transcripts in file
   >>> reader = BED_Reader(open("merlin_orfs.bed"),return_type=Transcript)

   >>> # do something with transcripts. here we just look at their names & attribute dictionaries
   >>> for transcript in reader:
   >>>     print(transcript.get_name() + ":\t" + str(transcript.attr))
   ORFL1W_(RL1):	{'cds_genome_end': 2299, 'color': '#000000', 'score': 0.0, 'cds_genome_start': 1366, 'type': 'mRNA', 'ID': 'ORFL1W_(RL1)'}
   ORFL2C:	{'cds_genome_end': 2723, 'color': '#000000', 'score': 0.0, 'cds_genome_start': 2501, 'type': 'mRNA', 'ID': 'ORFL2C'}
   ORFL3C:	{'cds_genome_end': 3015, 'color': '#000000', 'score': 0.0, 'cds_genome_start': 2934, 'type': 'mRNA', 'ID': 'ORFL3C'}
   [rest of output omitted]


|SegmentChains| and |Transcripts| can convert coordinates between the transcript
and the genome::

   >>> # load transcripts into a dictionary keyed on transcript ID
   >>> transcript_dict = { X.get_name() : X for X in BED_Reader(open("merlin_orfs.bed"),return_type=Transcript) }

   >>> # we'll use the two-exon, minus-strand gene ORFL83C as an example
   >>> demo_tx = transcript_dict["ORFL83C_(UL29)"]
   >>> demo_tx
   <Transcript segments=2 bounds=merlin:35004-37402(-) name=ORFL83C_(UL29)>

   >>> # get genomic coordinate of 1124th nucleotide from 5' end of ORFL83C
   >>> # right before the splice junction
   >>> demo_tx.get_genomic_coordinate(1124)
   ('merlin', 36277, '-')
    
   >>> # get genomic coordinate of 1125th nucleotide from 5' end of ORFL83C
   >>> # right after the splice junction
   >>> demo_tx.get_genomic_coordinate(1125)
   ('merlin', 36130, '-')

   >>> # and the inverse operation also works
   >>> demo_tx.get_segmentchain_coordinate("merlin",36130,"-")
   1126

.. _tour-get-counts:

|SegmentChains| can fetch vectors of data covering each position in the chain
from the 5' to 3' end (relative to the chain) from |GenomeArrays| (themselves
explained :ref:`below <tour-genome-array>`). For example, to count how many 5'
ends of sequencing reads appear at each position in a chain::

   >>> from plastid import BAMGenomeArray, FivePrimeMapFactory

   >>> # load read alignments, and map them to 5' ends
   >>> alignments = BAMGenomeArray(["SRR609197_riboprofile.bam"])
   >>> alignments.set_mapping(FivePrimeMapFactory())

   >>> # fetch the number of 5' ends of alignments at positions 300-320
   >>> demo_tx.get_counts(alignments)[320:340]
   array([  23.,    3.,   17.,   67.,   22.,    5.,   15.,   14.,   99.,
            26.,   13.,   27.,  112.,   34.,    1.,   13.,    0.,    4.,
             2.,   11.])

It is also possible to fetch sub-sections of a |Transcripts| or |SegmentChains|
as new |SegmentChains|::

   >>> # take first 200 nucleotides of  mRNA
   >>> subchain = demo_tx.get_subchain(0,200)
   >>> subchain
   <SegmentChain intervals=1 bounds=merlin:37202-37402(-) name=ORFL83C_(UL29)_subchain>

|Transcript| includes several convenience methods to fetch 5' UTRs, coding regions,
and 3'UTRs from coding transcripts::

   >>> demo_tx.get_utr5()
   <SegmentChain intervals=1 bounds=merlin:37353-37402(-) name=ORFL83C_(UL29)_5UTR>

   >>> demo_cds = demo_tx.get_cds()
   >>> demo_cds
   <Transcript intervals=2 bounds=merlin:35104-37353(-) name=ORFL83C_(UL29)_CDS>


|SegmentChain| and its subclasses can also fetch their sequences from dictionaries
of strings or :class:`Bio.SeqRecord.SeqRecord` objects. These sequences will
automatically be spliced and reverse-complemented, as necessary::

   >>> from Bio import SeqIO
   >>> genome = SeqIO.to_dict(SeqIO.parse(open("merlin_NC006273-2.fa"),"fasta"))
   >>> demo_tx.get_cds().get_sequence(genome)
   'ATGTCCGGCCGTCGCAAGGGCTGCTCGGCGGCCACGGCGTCTTCCTCCTCGTCGTCGCCGCCGTCGCGCCTTCCTCTGCCTGGGCACGCGCGTCGGCCGCGTCGCAAACGCTGCTTGGTACCCGAGG...'  
   # rest of output omitted


|SegmentChains| and |Transcripts| can do a lot more. For complete documentation
see |SegmentChain| and |Transcript| in :py:mod:`plastid.genomics.roitools`.
    
-------------------------------------------------------------------------------

.. _tour-genome-array:

|GenomeArray| & its subclasses
..............................
|GenomeArrays| are dictionary-like objects that  map quantitative data,
:term:`counts`, or :term:`read alignments`, to genomic positions.
Data can be imported from count files (`Wiggle`_, `bedGraph`_)
or alignment files (`bowtie`_ or `BAM`_ formats). For very large genomes a
sparse implementation is provided by |SparseGenomeArray|. A |BAMGenomeArray|
is provided for :term:`read alignments` in `BAM`_ format.

|GenomeArrays| can be indexed by |GenomicSegments| or |SegmentChains|. 
Doing so returns a vector of counts at each position in the |GenomicSegment|
or |SegmentChain|, with 5' to 3' coordinates relative to the chain (i.e.
for reverse-strand features, position 0 of the vector corresponds to
`segment.end`)::

   >>> # genomic segment
   >>> seg = GenomicSegment("merlin",1500,1600,"+")
   >>> alignments[seg]
   array([ 0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  1.,  1.,
       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
       4.,  1.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
       0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
       0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,
       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.,  0.,
       1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  0.,
       0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

   >>> # segment chain
   >>> alignments[demo_cds][:100]
   array([ 24.,   4.,   0.,   1.,   6.,   1.,   0.,   1.,  16.,   2.,   1.,
        1.,   2.,  13.,  17.,  13.,  13.,   2.,   3.,  23.,  10.,  39.,
       22.,  23.,  31.,  34.,  11.,  20.,  15.,   2.,   8.,  10.,   4.,
       11.,   9.,   5.,   5.,   4.,  13.,   5.,   2.,   0.,   2.,   4.,
        0.,   7.,  48.,  10.,  14.,   2.,   2.,   4.,   3.,   8.,   9.,
        0.,   9.,   8.,   8.,   9.,  10.,   9.,  14.,   3.,   9.,  33.,
        3.,   6.,  38.,   7.,   1.,  14.,   3.,  32.,  55.,  11.,   1.,
        4.,   1.,   9.,   9.,   1.,   3.,   2.,   0.,   6.,  17.,  21.,
        1.,  32.,   6.,   3.,  11.,   3.,   2.,   7.,  10.,   0.,  36.,
        4.])

   >>> # has same effects as calling the 'get_counts()' method
   >>> (demo_cds.get_counts(alignments) == alignments[demo_cds]).all()
   True

When importing :term:`read alignments`, users can specify a :term:`mapping rule`
to determine the genomic position(s) at which each alignment should be counted.
:data:`plastid` already includes mapping functions to map :term:`read alignments`:

  - to their fiveprime or threeprime ends, with or without offsets from
    that end (e.g. for :term:`P-site mapping <P-site offset>` for
    :term:`ribosome profiling data`)
     
  - fractionally over their entire lengths (for visualizing the full extent
    of transcripts in :term:`RNA-seq` data)
   
  - fractionally to all positions covered by a central portion of the read
    alignment, after excluding a user-defined number of positions on each
    send of the read (as in ribosome profiling data from *E. coli*
    :cite:`Oh2011` or *D. melanogaster* :cite:`Dunn2013`)

:term:`mapping rules <mapping function>` for |BAMGenomeArrays| can be changed
at runtime::

   >>> from plastid import FivePrimeMapFactory, ThreePrimeMapFactory
    
   >>> alignments.set_mapping(FivePrimeMapFactory())
   >>> demo_tx.get_cds().get_counts(alignments)[:50]
   array([ 24.,   4.,   0.,   1.,   6.,   1.,   0.,   1.,  16.,   2.,   1.,
            1.,   2.,  13.,  17.,  13.,  13.,   2.,   3.,  23.,  10.,  39.,
           22.,  23.,  31.,  34.,  11.,  20.,  15.,   2.,   8.,  10.,   4.,
           11.,   9.,   5.,   5.,   4.,  13.,   5.,   2.,   0.,   2.,   4.,
            0.,   7.,  48.,  10.,  14.,   2.])

   >>> # change to mapping with 15 nucleotide offset from 5' end
   >>> alignments.set_mapping(FivePrimeMapFactory(offset=15))
   >>> demo_tx.get_cds().get_counts(alignments)[:50]
   array([  3.,  26.,   8.,  17.,   4.,   4.,   9.,  27.,   7.,   3.,  17.,
           10.,  18.,  20.,   1.,  24.,   4.,   0.,   1.,   6.,   1.,   0.,
            1.,  16.,   2.,   1.,   1.,   2.,  13.,  17.,  13.,  13.,   2.,
            3.,  23.,  10.,  39.,  22.,  23.,  31.,  34.,  11.,  20.,  15.,
            2.,   8.,  10.,   4.,  11.,   9.])

   >>> # change to mapping from 3' end, with no offset
   >>> alignments.set_mapping(ThreePrimeMapFactory())
   >>> demo_tx.get_cds().get_counts(alignments)[:50]
   array([  5.,   6.,  14.,  17.,  24.,   5.,  14.,  19.,   4.,   5.,  11.,
            6.,   4.,   4.,   0.,   2.,   0.,   0.,   1.,   6.,  14.,  26.,
           25.,   4.,  23.,   7.,   8.,  24.,  11.,  11.,  22.,   9.,  14.,
            2.,   0.,   1.,   5.,   9.,   7.,   1.,   6.,   3.,   1.,   4.,
            5.,  15.,  15.,   6.,  17.,   8.])



|GenomeArrays| and subclasses can be exported to `wiggle`_ or `bedGraph`_
files for use in a :term:`genome browser`::

   >>> # export minus strand as a bedgraph file
   >>> with open("alignments_rc.wig","w") as fout:
   >>>     alignments.to_bedgraph(fout,"my_trackname","-")


`wiggle`_ or `bedGraph`_ files can be also imported into a |GenomeArray|
using the :meth:`~plastid.genomics.genome_array.GenomeArray.add_from_wiggle`
method::

   >>> new_data = GenomeArray()
   >>> new_data.add_from_wiggle(open("alignments_rc.wig"),"-")
    
   >>> demo_tx.get_cds().get_counts(new_data)[:50]
   array([  5.,   6.,  14.,  17.,  24.,   5.,  14.,  19.,   4.,   5.,  11.,
            6.,   4.,   4.,   0.,   0.,   0.,   0.,   0.,   6.,  14.,  26.,
           25.,   4.,  23.,   7.,   8.,  24.,  11.,  11.,  22.,   9.,  14.,
            2.,   0.,   0.,   5.,   9.,   7.,   1.,   6.,   3.,   1.,   4.,
            5.,  15.,  15.,   6.,  17.,   8.])

For further information, see:

 - The module documentation for :py:mod:`~plastid.genomics.genome_array`

 - In-depth discussion of :doc:`mapping rules <concepts/mapping_rules>`

-------------------------------------------------------------------------------

.. _tour-genome-hash:

|GenomeHash|, |BigBedGenomeHash|, and |TabixGenomeHash|
.......................................................

Often one needs to know whether any features overlap a specific region in the
genome, for example, to find transcripts that overlap one another.

But, it would be inefficient to scan an entire file to find overlapping features,
or to test whether two features overlap if we already know from their genomic
coordinates that they cannot.


|GenomeHash| and its subclasses avoid this problem by indexing features
by location. A |GenomeHash| may be created from a list or dictionary of features
(e.g. |SegmentChains| or |Transcripts|) in memory, or directly loaded from a
genome annotation (in `BED`_, `GTF2`_, `GFF3`_, or `PSL`_ format)::

   >>> from plastid import GenomeHash 
   >>> my_hash = GenomeHash(transcript_dict)
 
Having made a |GenomeHash|, we can ask what is where in the genome. For
example, to find all features between bases 10000-20000 on the plus
strand of chromosome *chrI*::

   >>> roi = GenomicSegment("merlin",10000,20000,"+")
   >>> my_hash[roi]
   [<Transcript segments=1 bounds=merlin:14307-14957(+) name=ORFL35W_(UL5)>,
    <Transcript segments=1 bounds=merlin:16522-17040(+) name=ORFL40W_(UL8)>,
    <Transcript segments=1 bounds=merlin:15814-16632(+) name=ORFL37W_(UL7)>,
    <Transcript segments=1 bounds=merlin:19793-21178(+) name=ORFL46W.iORF2>,
    <Transcript segments=1 bounds=merlin:12684-12929(+) name=ORFL25W>,
    <Transcript segments=1 bounds=merlin:13185-13406(+) name=ORFL30W>,
    <Transcript segments=1 bounds=merlin:19559-21178(+) name=ORFL46W>,
    <Transcript segments=1 bounds=merlin:9799-11193(+) name=ORFL23W_(RL12)>,
    <Transcript segments=1 bounds=merlin:13561-13779(+) name=ORFL33W>,
    <Transcript segments=1 bounds=merlin:12872-13192(+) name=ORFL26W>,
    <Transcript segments=1 bounds=merlin:18591-19559(+) name=ORFL45W_(UL11)>,
    <Transcript segments=1 bounds=merlin:19607-21178(+) name=ORFL46W.iORF1_(UL13)>,
    <Transcript segments=1 bounds=merlin:19053-19559(+) name=ORFL45W.iORF1>,
    <Transcript segments=1 bounds=merlin:18467-18685(+) name=ORFL44W>,
    <Transcript segments=1 bounds=merlin:17867-18142(+) name=ORFL42W>,
    <Transcript segments=1 bounds=merlin:14914-15906(+) name=ORFL36W_(UL6)>,
    <Transcript segments=1 bounds=merlin:13770-14369(+) name=ORFL34W_(UL4)>,
    <Transcript segments=1 bounds=merlin:11138-12169(+) name=ORFL24W_(RL13)>,
    <Transcript segments=1 bounds=merlin:14565-14957(+) name=ORFL35W.iORF1>]

Or on both strands::

   >>> my_hash.get_overlapping_features(roi,stranded=False)
   [<Transcript segments=1 bounds=merlin:19793-21178(+) name=ORFL46W.iORF2>,
    <Transcript segments=1 bounds=merlin:17609-17857(-) name=ORFL43C.iORF1>,
    <Transcript segments=1 bounds=merlin:19607-21178(+) name=ORFL46W.iORF1_(UL13)>,
    <Transcript segments=1 bounds=merlin:17609-17968(-) name=ORFL43C>,
    <Transcript segments=2 bounds=merlin:7811-13052(-) name=ORFL27C>,
    <Transcript segments=1 bounds=merlin:15643-15915(-) name=ORFL38C>,
    <Transcript segments=1 bounds=merlin:13561-13779(+) name=ORFL33W>,
    <Transcript segments=1 bounds=merlin:19053-19559(+) name=ORFL45W.iORF1>,
    <Transcript segments=1 bounds=merlin:17867-18142(+) name=ORFL42W>,
    <Transcript segments=1 bounds=merlin:14914-15906(+) name=ORFL36W_(UL6)>,
    <Transcript segments=1 bounds=merlin:14307-14957(+) name=ORFL35W_(UL5)>,
    <Transcript segments=1 bounds=merlin:16522-17040(+) name=ORFL40W_(UL8)>,
    <Transcript segments=1 bounds=merlin:15698-15937(-) name=ORFL39C>,
    <Transcript segments=1 bounds=merlin:17061-17294(-) name=ORFL41C>,
    <Transcript segments=1 bounds=merlin:13169-13402(-) name=ORFL31C.iORF1>,
    <Transcript segments=1 bounds=merlin:13185-13406(+) name=ORFL30W>,
    <Transcript segments=1 bounds=merlin:19559-21178(+) name=ORFL46W>,
    <Transcript segments=1 bounds=merlin:12872-13192(+) name=ORFL26W>,
    <Transcript segments=1 bounds=merlin:12957-13226(-) name=ORFL29C>,
    <Transcript segments=1 bounds=merlin:13770-14369(+) name=ORFL34W_(UL4)>,
    <Transcript segments=1 bounds=merlin:12910-13140(-) name=ORFL28C.iORF1>,
    <Transcript segments=1 bounds=merlin:14565-14957(+) name=ORFL35W.iORF1>,
    <Transcript segments=1 bounds=merlin:15814-16632(+) name=ORFL37W_(UL7)>,
    <Transcript segments=1 bounds=merlin:12684-12929(+) name=ORFL25W>,
    <Transcript segments=1 bounds=merlin:9799-11193(+) name=ORFL23W_(RL12)>,
    <Transcript segments=1 bounds=merlin:13110-13439(-) name=ORFL31C_(UL2)>,
    <Transcript segments=1 bounds=merlin:18591-19559(+) name=ORFL45W_(UL11)>,
    <Transcript segments=1 bounds=merlin:11138-12169(+) name=ORFL24W_(RL13)>,
    <Transcript segments=1 bounds=merlin:12899-13159(-) name=ORFL28C>,
    <Transcript segments=1 bounds=merlin:18467-18685(+) name=ORFL44W>]

Does anything interesting overlap *ORFL83C_(UL29)*?

.. code-block:: python

   >>> my_hash[demo_tx]
   [<Transcript segments=1 bounds=merlin:37077-37449(-) name=ORFL84C>,
    <Transcript segments=1 bounds=merlin:33081-35058(-) name=ORFL79C_(UL27)>,
    <Transcript segments=1 bounds=merlin:37382-37898(-) name=ORFL85C_(UL30)>,
    <Transcript segments=2 bounds=merlin:35004-37403(-) name=ORFL83C_(UL29)>]

For more information, see the module documentation for :mod:`~plastid.genomics.genome_hash`.

-------------------------------------------------------------------------------


See also
--------
  - :doc:`examples` 

  - Detailed :ref:`module documentation <modindex>` for complete descriptions
    of the attributes and methods of these and other data structures

      - :mod:`plastid.genomics.roitools`

      - :mod:`plastid.genomics.genome_array`

      - :mod:`plastid.genomics.genome_hash`

      - :mod:`plastid.readers`
