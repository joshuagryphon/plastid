Manipulating genomic sequence
=============================
In this tutorial we will fetch the sequences of genomic `features <feature>`_,
using the :doc:`/test_dataset`.


Fetching sequence from |SegmentChains|
--------------------------------------
|SegmentChains| and |Transcripts| can fetch their own genomic sequence via
their :meth:`~yeti.genomics.roitools.SegmentChain.get_sequence` method.
As an argument, the method expects a dictionary-like object of string-like
objects, where the dictionary keys correspond to chromosome names, and
the strings correspond to the forward-strand sequence.

    >>> little_genome = { "chrA" : "NNNNNNNNNNTCTAGACGATACANNNNNNNNNNCTACGATA" }

Sequences are fully-spliced, if the |SegmentChain| contains one or more
|GenomicSegment|::

    >>> from yeti.genomics.roitools import GenomicSegment, SegmentChain

    >>> little_chain = SegmentChain(GenomicSegment("chrA",10,23,"+"),
    >>>                             GenomicSegment("chrA",33,41,"+"),
    >>>                             ID="little_chain")
    >>> little_chain.get_sequence(little_genome)
    'TCTAGACGATACACTACGATA'
    

For reverse-strand features, sequences are automatically reverse-complemented::

    >>> antisense_chain = little_chain.get_antisense()
    >>> antisense_chain
    <SegmentChain segments=2 bounds=chrA:10-41(-) name=chrA:10-23^33-41(-)>

    >>> antisense_chain.get_sequence(little_genome)
    'TATCGTAGTGTATCGTCTAGA'


For convenience, |SegmentChain| and |Transcript| also implement
:meth:`~yeti.genomics.roitools.SegmentChain.get_fasta`, which formats output
in `FASTA`_ format::

    >>> print(little_chain.get_fasta(little_genome))
    >little_chain
    TCTAGACGATACACTACGATA


Supported sequence file formats
-------------------------------
Genomic sequence can be stored in a number of formats. Two of the most common are:

  - `FASTA`_: a non-indexed, text-base file of one or more sequences.

  - `2bit <twobit>`_: an indexed, binary format created for large (e.g. mammalian)
     genomes. Beacause `2bit <twobit>`_ files are binary and indexed, they use
     far less memory than `FASTA`_ files.

:data:`yeti` doesn't implement any of its own file readers, but it is compatible
with any parser that returns sequences as a dictionary-like object of string-like
objects. This allows users to use implementations in `Biopython`_ (for `FASTA`_, 
EMBL, & GenBank formats) and `twobitreader`_ (for `2bit <twobit>`_) files.
For example, to load a FASTA file using `Biopython`_::

    >>> # load CMV genome from test dataset
    >>> from Bio import SeqIO
    >>> genome = SeqIO.to_dict(SeqIO.parse(open("merlin_NC006273-2.fa"),"fasta"))
    >>> genome
    {'merlin': SeqRecord(seq=Seq('CCATTCCGGGCCGTGTGCTGGGTCCCCGAGGGGCGGGGGGGTGTTTTCTGCGGG...GCT', SingleLetterAlphabet()), id='merlin', name='merlin', description='merlin gi|155573622|ref|NC_006273.2| Human herpesvirus 5 strain Merlin, complete genome', dbxrefs=[])}

In `Biopython`_, nucleic acids are represented as `SeqRecord`_ objects
rather than strings. |SegmentChain| and |Transcript| don't mind::

    >>> # load CMV annotations from test dataset
    >>> from yeti.genomics.roitools import Transcript
    >>> from yeti.readers.bed import BED_Reader
    >>> transcripts = list(BED_Reader(open("merlin_orfs.bed"),return_type=Transcript))
    >>> transcripts[0]
    <Transcript segments=1 bounds=merlin:1316-2398(+) name=ORFL1W_(RL1)>

    >>> transcripts[0].get_sequence(genome)[:200]
    'GCTCGCCTATTTAACCTCCACCCACTACAACACACACATGCCGCACAATCATGCCAGCCACAGACACAAACAGCACCCACACCACGCCGCTTCACCCAGAGGACCAACACACGTTACCCTTACACCACAGCACCACACAACCTCATGTCCAAACTTCGGACAAACACGCCGACAAACAACACCGCACGCAGATGGAGCTC'


Similarly, `TwoBitFile`_ objects from `twobitreader`_  can be directly passed
to :meth:`~yeti.genomics.roitools.SegmentChain.get_sequence`, because they 
inherit from :class:`dict` and return strings::

# load CMV genome as a 2bit file
    >>> from twobitreader import TwoBitFile
    >>> twobit_genome = TwoBitFile("merlin_NC006273-2.2bit")
    >>> twobit_genome.keys()
        ['merlin']

    >>> transcripts[0].get_sequence(twobit_genome)[:200]
    'GCTCGCCTATTTAACCTCCACCCACTACAACACACACATGCCGCACAATCATGCCAGCCACAGACACAAACAGCACCCACACCACGCCGCTTCACCCAGAGGACCAACACACGTTACCCTTACACCACAGCACCACACAACCTCATGTCCAAACTTCGGACAAACACGCCGACAAACAACACCGCACGCAGATGGAGCTC'


Manipulating sequence
---------------------
Tools for further manipulating sequence (e.g. reverse-complementing, translating)
are supplied in `Biopython`_'s `Seq`_ and `SeqRecord`_ objects::

    >>> # SeqRecord examples
    >>> from Bio.Alphabet import generic_dna
    >>> from Bio.Seq import Seq

    >>> seq = Seq(transcripts[0].get_cds().get_sequence(genome),generic_dna)
    >>> seq.translate()
    Seq('MPATDTNSTHTTPLHPEDQHTLPLHHSTTQPHVQTSDKHADKQHRTQMELDAAD...PW*', HasStopCodon(ExtendedIUPACProtein(), '*'))

Fuller explanations and further examples can be found in the `Biopython`_
documentation for `Seq`_ and `SeqRecord`_.

-------------------------------------------------------------------------------

See also
--------
  - `Biopython`_ documentation for manipulation of nucleic acid sequences
  - `twobitreader`_ documentation
  - `UCSC file format FAQ`_ for details on sequence file formats