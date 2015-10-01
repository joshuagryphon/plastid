Related resources
=================

If you are interested in :data:`plastid`, the following projects might also
be helpful for you:


Similar projects, for interactive analysis of sequencing data
-------------------------------------------------------------
`HTSeq`_
    A Python package designed for analysis of high-throughput
    sequencing data.

`metaseq`_
    A Python package for analysis of genomics data. It is very
    similar in intent to :data:`plastid` in that it introduces
    simple, unified APIs to access genomic data of many file
    types. In addition, it includes nice plotting utilities
    for interactive analysis

`riboseqr`_
    A toolkit for analysis of :term:`ribosome profiling` data,
    written in `R`_. It implements many standard workflows.
   
`ribogalaxy`_
    A web-based platform for analysis of :term:`ribosome profiling`
    data, integrating `riboseqr`_, `galaxy`_, and other tools.


Gene expression analysis
------------------------
`DESeq`_
    Statistical models for assessment of differential gene expression,
    applicable to RNA-seq, :term:`ribosome profiling`, and many other
    types of :term:`high-throughput sequencing` data
    
    `DESeq`_ can be used to test for significant differences in expression
    counts obtained using the :mod:`~plastid.bin.cs` or
    :mod:`~plastid.bin.counts_in_region` scripts

`cufflinks`_
    A software suite for transcript assembly and differential expression
    analysis of RNA-seq data
 
`kallisto`_
    Software for measurement of gene expression from RNA-seq data


General-purpose manipulation
----------------------------
`Samtools`_
    Manipulate read alignments in SAM and `BAM`_ files on the command line
 
`Pysam`_
    Python bindings for `Samtools`_. :data:`plastid` uses `Pysam`_ internally
    for parsing of `BAM`_ and `tabix`_-compressed files

`bedtools`_
    Fast command-line tools that perform arithmetic on annotations of continuous
    genomic features, and that count read coverage and/or other properties
    for those regions

`pybedtools`_
    Python bindings for `bedtools`_

`Jim Kent's utilities`_
    Convert text-based genomic files to randomly accessible, indexed binary 
    formats (e.g. `BED`_ to `BigBed`_, `wiggle`_ and `bedGraph`_
    to `BigWig`_, `FASTA`_ to `2bit <twobit>`_, et c)


Genome browsers
---------------
`Integrative Genome Viewer <IGV>`_
    A lightweight and versatile genome browser created
    by the `Broad Institute <www.broadinstitute.org>`_. `IGV`_ is suitable
    for laptops & desktops.

`UCSC Genome Browser`_
    A web-based genome browser developed by University of California,
    Santa Cruz. The `UCSC Genome Browser`_ integrates with UCSC's large
    database of genomes, annotations, and tracks of quantitive data.
    It also offers many tools for visualization and manipulation
    of genomics data.
