 .. glossary ::
    :sorted:

    alignment
    read alignment
        A record matching a sequencing read to the genomic coordinates from
        which it presumably derived. These are produced by running sequencing
        data through alignment programs, such as `Bowtie`_, `Tophat`_, or `BWA`_.
        The most common format for short read alignments is `BAM`_.

    annotation
        A file that describes locations of features (e.g. genes, mRNAs)
        in a genome. These come in various formats, e.g.  `BED`_, `BigBed`,
        `GTF2`_, `GFF3`_, `PSL`_.

    count file
        A file that assigns quantitative data -- for example, read alignment
        counts, or conservation scores -- to genomic coordinates. Strictly
        speaking, these include  `bedGraph`_ or `wiggle`_ files but ``yeti``
        can also treat :term:`alignment` files in `Bowtie`_ or `BAM`_ format
        as count files, if a :term:`mapping rule` is applied.

    crossmap
        A genomic annotation that identifies regions of the genome that
        cannot give rise to uniquely-mapping reads due to repetitive sequence.
        Crossmaps are functions of genome sequence, read length, and alignment
        parameters (e.g. the number of mismatches allowed). These may be
        generated from genome sequence using the included
        :py:mod:`~yeti.bin.crossmap` script.

    feature
        A region of the genome with interesting or specific properties, such
        as a gene, an mRNA, an exon, a centromere, et c.

    mapping rule
    read mapping rule
        A function that describes how a read alignment is to be mapped
        to the genome for positional analyses. Reads typically are mapped
        to their fiveprime or threeprime ends, with an offset of 0 or more
        nucleotides that can optionally depend on the read length.
        
        For example, to map a ribosome-protected mRNA fragment to the nucleotide
        position corresponding to the A-site of the ribosome that protected it,
        (i.e. a nucleotide within the codon being decoded), it is common to use
        a 12 nucleotide offset from the 3' end of the fragment.

    metagene
    metagene average
        An average of quantitative information over genes aligned at some
        internal feature. For example, an average of ribosome densit across
        all genes, aligned at their start codons. Or, perhaps, an average
        across all genes of nucleotide sequence conservation across the 12
        fly genomes surrounding 5' splice sites of first introns. See the
        documentation for the `:py:mod:~yeti.bin.metagene` script for more
        explanation.

    footprint
    ribosome-protected fragment
        A fragment of mRNA protected from nuclease digestion by a ribosome
        during ribosome profiling or other molecular biology assays.

    roi
    region of interest
        A region of the genome or of a transcript that contains an interesting
        :term:`feature`.

    RPKM
        Reads per kilobase per million reads in a dataset. This is a unit of
        sequencing density that is normalized by sequencing depth (in millions of
        reads) and by the length of the region of interest (in kb).

    single-end sequencing
        A high-throughput sequencing technique that generates short reads
        of approximately 50-100 nt in length

    paired-end sequencing
        A high-throughput sequencing technique in which 50-100 nucleotides
        of each end of a ~300 nucleotide sequence are read, and reported
        as a pair.

