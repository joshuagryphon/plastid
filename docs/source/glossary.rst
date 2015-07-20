Glossary of terms
=================

 .. glossary ::
    :sorted:

    alignment
    read alignments
        A record matching a short sequence of DNA or RNA to a region of identical or similar
        sequence in a genome. In a :term:`high-throughput sequencing` experiment,
        alignment of short reads identifies the genomic coordinates from which
        each read presumably derived.
        
        These are produced by running sequencing
        data through alignment programs, such as `Bowtie`_, `Tophat`_, or `BWA`_.
        The most common format for short read alignments is `BAM`_.

    annotation
        A file that describes locations and properties of :term:`features <feature>`
        (e.g. genes, mRNAs, SNPs, start codons) in a genome. Annotation files
        come in various formats, such as `BED`_, `BigBed`_, `GTF2`_, `GFF3`_,
        and `PSL`_, among others. In a :term:`high-throughput sequencing`
        experiment, it is essential to make sure that the coordinates in the
        :term:`annotation` correspond to the :term:`genome build` used
        to generate the alignments.

    counts
        Colloquially, the number of :term:`read alignments` overlapping a region
        of interest, or mapped to a nucleotide.
    
    count file
        A file that assigns quantitative data -- for example, read alignment
        counts, or conservation scores -- to genomic coordinates. Strictly
        speaking, these include  `bedGraph`_ or `wiggle`_ files but :py:obj:`yeti`
        can also treat :term:`alignment` files in `bowtie`_ or `BAM`_ format
        as count files, if a :term:`mapping rule` is applied.

    crossmap
    mask file
    mask annotation file
        A genomic :term:`annotation` that identifies regions of the genome that
        cannot give rise to uniquely-mapping reads due to repetitive sequence.
        :term:`mask files <mask file>` are functions of genome sequence, read length, and alignment
        parameters (e.g. the number of mismatches allowed). These may be
        generated from genome sequence using the included
        :py:mod:`~yeti.bin.crossmap` script.

    factory function
        A function that produces functions

    feature
        A region of the genome with interesting or specific properties, such
        as a gene, an mRNA, an exon, a centromere, et c.

    genome assembly
    genome build
        A specific edition of a genome sequence for a given organism. These
        are updated over time as sequence data is added and/or corrected.
        When an assembly is updated, frequently the lengths of the chromosomes or
        contigs change as sequences are corrected. 

    genome browser
        Software used for visualizing genomic sequence, :term:`feature`
        annotations, :term:`read alignments`, and other quantitative data
        (e.g. nucleotide-wise sequence conservation). Popular genome browsers
        include `IGV`_ and the `UCSC genome browser`_. 

    deep sequencing
    high-throughput sequencing
        A group of experimental techniques that produce as output millions of
        reads (strings) of short DNA sequences.

    k-mer
        A sequence *k* nucleotides long.

    mapping rule
    mapping function
        A function that describes how a read alignment is mapped
        to the genome for positional analyses. Reads typically are mapped
        to their fiveprime or threeprime ends, with an offset of 0 or more
        nucleotides that can optionally depend on the read length.
        
        For example, ribosome-protected mRNA fragments are frequently mapped
        to their :term:`P-site offset` by using a 15 nucleotide offset 
        from the threeprime end of the fragment.

        See :doc:`/concepts/mapping_rules` for more information.

    metagene
    metagene average
        An average of quantitative data over one or more
        genomic regions (often genes or transcripts) aligned at some internal feature.
        For example, a :term:`metagene` profile could be built around:
      
          - the average of ribosome density surrounding the start codons of all 
            transcripts in a :term:`ribosome profiling` dataset
        
          - an average phylogenetic conservation score surounding the 5' splice
            site of the first introns of all transcripts
      
        See :doc:`/examples/metagene` and/or the module documentation for the
        :py:mod:`~yeti.bin.metagene` script for more explanation.

    footprint
    ribosome-protected footprint
        A fragment of mRNA protected from nuclease digestion by a ribosome
        during ribosome profiling or other molecular biology assays.

    ribosome profiling
        A :term:`high-throughput sequencing` technique that captures the positions
        of all ribosomes on all RNAs at a snapshot in time. See :cite:`Ingolia2009`
        for more details

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
        of approximately 50-100 nt in length.

    paired-end sequencing
        A high-throughput sequencing technique in which 50-100 nucleotides
        of each end of a ~300 nucleotide sequence are read, and reported
        as a pair.

    P-site offset
        Distance from the 5' or 3' end of a ribosome-protected footprint
        to the P-site of the ribosome that generated the footprint (see
        :cite:`Ingolia2009`, fig. 2B). Because the P-site is the site where
        peptidyl elongation occurs, read alignments from :term:`ribosome profiling`
        are frequently mappped to their P-site offsets, as opposed to their 5'
        or 3' ends.
        
        P-site offsets may be estimated from ribosome profiling data
        using the :py:mod:`~yeti.bin.psite` script.

    start codon peak
        Large peaks of :term:`ribosome-protected footprint` visible over initiator codons
        in ribosome profiling data. These arise because the kinetics of
        translation initiation are slow compared to the kinetics of
        elongation, causing a build-up over the initiator codon.

    stop codon peak
        Large peaks of :term:`ribosome-protected footprint` visible
        over stop codons in some ribosome profiling datasets. These
        arise because the kinetics of translation termination are 
        slow compared to the kinetics of elongation, causing a build-up
        over termination codons. These peaks are frequently absent
        from datasets if tissues are pre-treated with elongation
        inhibitors (e.g. cycloheximide) before lysis and sample prep.

    sub-codon phasing
    triplet periodicity
        A feature of :term:`ribosome profiling` data. Because ribosomes
        step three nucleotides in each cycle of translation elongation,
        in many :term:`ribosome profiling` datasets a triplet periodicity
        is observable in the distribution of
        :term:`ribosome-protected footprints <footprint>`, in which 70-90%
        of the reads on a codon fall within the first of the three codon
        positions. This allows deduction of translation reading frames,
        if the reading frame is not known *a priori.* See :cite:`Ingolia2009`
        for more details
