A simple alignment and quantitation workflow
============================================

This tutorial covers a rudimentary workflow for aligning reads and doing some
QC on the dataset. It assumes that the genome has been set up, as described in 
:doc:`/examples/a1_genome_setup`.

.. note::

   This tutorial is very much under consruction.
   
This document includes the following sections:

.. contents::
   :local:


.. _starting-out-aligments:

Preflight checks
----------------

Remove cloning adaptors, if present
...................................

This stage is now frequently taken care of by sequencing facilities. But, if
not, you may need to remove untrimmed adaptor sequences manually.

Good options for this include: 

 - the ``fastx_clipper`` utility from the `fastx toolkit`_
   
 - `Trimmomatic`_, which has nice support for paired-end reads

Exact options will depend upon the cloning strategies used in your library prep.
Please refer to their documentation as needed.


Check raw read quality
......................

A number of packages exist for checking the quality of the raw data themselves,
and are useful, for example, for spotting a failed cycle of sequencing. 
The Illumina pipeline creates HTML reports for a number of useful metrics during
basecalling. Your sequencing facility can provide these to you.

If not, `FastQC`_ produces useful plots of base-wise quality scores, read length
distributions, and nucleotide distributions. Its documentation contains 
examples of what you should expect to see if things work out.



Perform alignments
------------------

Prealignment against rRNA
.........................

If working on a sample derived from RNA (RNA-seq, :term:`ribosome profiling`, 
et c), it is frequently useful to computationally filter out fragments derived
from rRNA. Because these are abundant -- even *after* rRNA depletion or polyA
selection -- removing them saves processing time and
disk space downstream.

`bowtie`_ is an excellent tool for filtering out rRNA-derived fragments. To use
it, one must first obtain or construct a `fasta`_ file of rRNA sequences from
the organism of interest. These may be obtained from `GenBank`_. Once
constructed (in this example, assume it is named `species_rRNA.fa`), build a
bowtie index (in this example, called `species_rRNA`). From the terminal:

.. code-block:: shell

   $ bowtie-build species_rRNA species_rRNA.fa

Then, align the data, discarding any reads mapping to the rRNA, and saving those
that do not align to a separate file, here called `my_reads_rRNA_unalign.fq`:

.. code-block:: shell 

   $ bowtie -p4 -v3 species_rRNA my_reads.fastq \
            --un my_reads_rRNA_unalign.fq >/dev/null

Reads that do not align (in `my_reads_rRNA_unalign.fq`), may then be used as 
input for additional stages of prealignment, or as input for chromosomal
alignment.


Other uses of (pre)alignment
............................

Depending upon how you intend to analyze your data, it may be useful to separate
out other classes of reads with additional stages of prealignment:

 - To remove expected sequence contaminants:
 
   For example, flies eat yeast. But, often, if we're performing an 
   experiment on flies, we aren't too concerned with which genes the 
   yeast are expressing. In this case, it is useful to make a `bowtie`_ index
   of the yeast transcriptome, and filter those reads out, just as we did 
   with the rRNA above.

 - To measure transposon expression:
 
   Because transposons exist throughout the genome in multiple, quite similar copies,
   they produce multimapping reads.
   
   If retrotransposon activity is of interest, one way to capture it is
   to make a database of consensus sequences of each type of transposon in your
   genome of interest, and to pre-align to it allowing a generous number of
   mismatches. Instead of discaring the reads (as we did with rRNA), retain the
   alignments, using them for expression quantitation. 

As previously, take all the reads which did not align in the pre-alignment
stages, and put them into chromosomal alignment, below.  



Chromosomal alignment
.....................

Having filtered out reads derived from rRNA, we can align the remaining data to
the genome. In this example, we use `Tophat`_, because it is capable of
performing spliced (and other gapped) alignments. We'll use splice junctions
from the ``juncs`` file made in :ref:`the genome setup tutorial <starting-out-juncs-file>`.

.. code-block:: shell

   # run tophat
   $ tophat -o my_chr_alignments \
            --bowtie1 \
            --raw-juncs splice_junctions.juncs \
            --no-novel-juncs \
            /path/to/chromosme/bowtie/index \
            my_reads_rRNA_unalign.fq 

   # rename output file
   $ mv my_chr_alignments/accepted_hits.bam chr_alignments.bam

   # index the BAM file for use in plastid, IGV, et c
   $ samtools index chr_alignments.bam


.. note::

   Selection of chromosomal alignment parameters is a complex topic, and beyond
   the scope of this tutorial.
   
   For more details, see manual for `Tophat`_ (or your favorite aligner)
   and other forums dedicated to this purpose. It is important to choose
   alignment parameters tailored to your own needs.



Analysis
--------


Special considerations for ribosome profiling
.............................................

If analyzing :term:`ribosome profiling` data, it is helpful to estimate
the location of the ribosomal P-site and to analyze the :term:`sub-codon phasing`
as quality control metrics.

To do so, make sure you have generated :term:`maximal spanning windows <maximal spanning window>` from
your genome annotation using the |metagene| script's ``generate`` subprogram
(see the :ref:`metagene tutorial <metagene-generate>` and the |metagene|
script documentation for details).

The output from ``metagene generate`` can then be used for P-site estimation
and analysis of :term:`sub-codon phasing`.


P-site estimation
"""""""""""""""""

Estimation of ribosomal :term:`P-site offsets <P-site offset>` is important for
all subsequent position-wise analyses in ribosomal profiling (discussed in 
:cite:`Ingolia2009`).

See :doc:`/examples/p_site` for background and step-by-step instructions.


Sub-codon phasing
"""""""""""""""""

After determining :term:`P-site offsets <P-site offset>`, it is possible to
examine the :term:`sub-codon phasing` found in your ribosome profiling data. See
:doc:`/examples/phasing` for step-by-step instructions and information on the
|phase_by_size| script.


Read counting, gene expression, and differential expression analysis
....................................................................

This is a large topic, and the details of how to do this depend upon your
experiment. One workflow we like is to:

 1. Use |cs| to create a corrected genome annotation
 
 2. Use |cs| to count the number of reads aligning to each gene in each dataset
 
 3. Combine relevant columns the output from |cs| from each sample into a single
    master table 

 4. Using the master table from (3), perform some clustering to make sure samples
    behave as expected. For example, replicates should cluster together, while
    differently treated samples  might not. We don't have a tutorial on this
    right now, but `DESeq2`_ has a fabulous tutorial on this on its home page,
    in its vignettes. We recommend you check it out!
     
 5. Use `DESeq2`_ to perform differential expression analysis and test for
    significance

A step-by-step discussion of a similar workflow, using the simpler
|counts_in_region| script instead of |cs| appears in  :doc:`/examples/gene_expression`.
At present, steps (1-3) and (5) are discussed.


Visualization in genome browsers
--------------------------------

Modern genome browsers can import `BAM`_ files of alignments directly. By
default, they tend to plot individual alignments and a summary track of read
coverage over each nucleotide.

Frequently it is useful to plot some aspect of the data, rather than raw read
coverage, as a function of nucleotide position in the genome. For example, in
:term:`ribosome profiling` data, it is useful to plot the total number of ribosomes
translating a particular nucleotide.

This data can be extracted from read alignments by use of :term:`mapping rules <mapping rule>`.
See :doc:`/concepts/mapping_rules` for a further discussion which
:term:`mapping rules <mapping rule>` are available, and how to use them.

To export transformed data as a browser track in `bedGraph`_ or `wiggle`_
formats, see the |make_wiggle| script. If working on a large (e.g. plant or
metazoan) genome, it might be helpful to convert the output from |make_wiggle|
into to a `BigWig`_ file, using `Jim Kent's utilities`_ from UCSC.

.. TODO: give example