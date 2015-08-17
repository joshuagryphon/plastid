Related resources
=================

If you are interested in :data:`yeti`, the following projects might also
be helpful for you:

 .. TODO finish this page

Similar projects, for interactive analysis of sequencing data
-------------------------------------------------------------
:data:`HTSeq`
    [description]

`metaseq`_
    [description]


Short read aligners
-------------------
`bowtie`_
    A fast short-read aligner for ungapped alignments.

`TopHat`_
    A short read aligner built atop `bowtie`_ that can perform gapped alignments
    and discover splice junctions.


Gene expression analysis
------------------------
`DESeq`_
    Statistical models for assessment of differential gene expression,
    applicable to RNA-seq, :term:`ribosome profiling`, and many other
    types of :term:`high-throughput sequencing` data.
    
    `DESeq`_ can be used to test for significant differences in expression
    counts obtained using the :mod:`~yeti.bin.cs` or
    :mod:`~yeti.bin.counts_in_region` scripts.

`cufflinks`_
    A software suite for transcript assembly and differential expression
    analysis of RNA-seq data.
 
`kallisto`_
    Software for measurement of gene expression from RNA-seq data.


General-purpose manipulation
----------------------------
`Samtools`_
    Manipulate read alignments in SAM and `BAM`_ files on the command line
 
`Pysam`_
    Python bindings for `Samtools`_. :data:`yeti` uses `Pysam`_ internally
    for parsing of `BAM`_ and `tabix`_-compressed files.

`bedtools`_
    Manipulate annotations of continuous (unspliced) genomic regions,
    and count read coverage for those regions

`Jim Kent's utilities`_
    Index and convert files (e.g. `BED`_ to `BigBed`_, `wiggle`_ and `bedGraph`_
    to `BigWig`_, `FASTA`_ to `2bit <twobit>`_, et c).


Genome browsers
---------------
`Integrative Genome Viewer <IGV>`_
    A lightweight and versatile genome browser created
    by the `Broad Institute <www.broadinstitute.org>`_. `IGV`_ is suitable
    for laptops & desktops.

`UCSC Genome Browser`_
    A web-based genome browser developed by University of California,
    Santa Cruz. The `UCSC Genome Browser`_ integrates with UCSC` large
    database of genomes and annotations, and offers many tools for manipulation
    or analysis of genomics data.
