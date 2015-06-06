Command-line scripts
====================
:py:data:`yeti` comes with a set of command-line scripts
that implement common workflows in RNA-seq and ribosome profiling.
See the individual titles below for in-depth descriptions and example usage:


Sequencing & ribosome profiling
-------------------------------

:py:mod:`~yeti.bin.crossmap`
	Empirically determine regions of the genome that give rise to repetitive sequence, 
	so that these may be excluded from downstream analysis

:py:mod:`~yeti.bin.cs3`
	Tabulate read densities (in RPKM and raw counts) over genes, 
	making several corrections to gene boundaries when genes overlap one another, 
	or overlap repetitive regions as determined by
	:py:mod:`~yeti.bin.crossmap` or other sources

:py:mod:`~yeti.bin.get_count_vectors`
	Fetch vectors of read counts at each position of each region of interest
	specified in an annotation file,  masking repetitive regions as determined by
	:py:mod:`~yeti.bin.crossmap` or other sources.

:py:mod:`~yeti.bin.make_wiggle`
	Create browser tracks in `Wiggle`_ or `bedGraph`_ format from alignments 
	in `bowtie`_ or `BAM`_ format.

:py:mod:`~yeti.bin.metagene`
	Create metagene profiles of sequencing data

:py:mod:`~yeti.bin.phase_by_size`
	Measure sub-codon phasing for ribosome profiling data, stratified by
	read length.
	
	**note:** output formats will change as this script is redesigned 

:py:mod:`~yeti.bin.psite`
	Estimate P-site offsets for ribosome profiling data various read lengths


Splice junctions
----------------

:py:mod:`~yeti.bin.findjuncs`
	Create a BED file of splice junctions in a transcript annotation

:py:mod:`~yeti.bin.make_tophat_juncs`
	Create a `tophat <http://ccb.jhu.edu/software/tophat/index.shtml>`_
	``.juncs`` file from a genome annotation

:py:mod:`~yeti.bin.mergejuncs`
	Create a non-redundant `BED`_ file of splice junctions from two or more
	annotations

:py:mod:`~yeti.bin.slidejuncs`
	Compare discovered splice junctions to known annotations, and classify them
	as truly novel, compatible with nearby known splice junctions, or buried
	in non-unique repetitive sequence


Miscellaneous
-------------
:py:mod:`~yeti.bin.gff_parent_types`
	Examine schema of parent-child relationships in a GFF3 file

	
 .. toctree::
    :maxdepth: 2
    :hidden: 
    :glob:
	
    generated/yeti.bin.*
