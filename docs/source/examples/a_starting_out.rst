Getting started with analysis
=============================

This tutorial is under construction. When complete, it will describe a sample
workflow for sequencing analysis, using :data:`plastid` and other tools.



.. contents::
   :local:


Download some tools
-------------------

In this tutorial, we use the following tools. They are useful to have on-hand
in many circumstances:

 - `samtools`_ to index and manipulate `BAM`_ files

 - `bowtie`_ or another short read aligner, and possibly also `TopHat`_

 - The Hannon Lab's `fastx toolkit`_ for removing cloning adaptors
 
 - `Jim Kent's utilities`_ for making `BigBed`_ or `BigWig`_ files



.. _starting-out-annotation:
 
Obtain a genome sequence & matching annotation
----------------------------------------------

TODO


.. _starting-out-derivative-files:

Create derivative files
-----------------------

TODO

.. It is often useful to pre-compute a number of files for use later:

   - A `BED`_/`BigBed`_ file of transcripts with gene IDs, rather than a `GTF2`_
     or `GFF3`_ annotation

   - A ``.juncs`` file for aligning data in `TopHat`_
   - A `bowtie`_ index
   - maximal spanning windows around start codon



.. _starting-out-aligments:

Perform alignments
------------------

Remove cloning adaptors, if present
...................................

If your pipeline does not automatically remove cloning adaptor sequences, you
may need to remove these manually. A good option for this is the
``fastx_clipper`` utility from the `fastx toolkit`_. See their documentation for
detailed instructions.



Prealignment against rRNA
.........................

If working on a sample derived from RNA, it is frequently useful to
computationally filter out fragments derived from rRNA. These tend to be present
in high quantities even in poly-A selected or rRNA subtracted samples. Filtering
these out will save processing speed and disk space downstream.

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

   $ bowtie -p4 -v3 species_rRNA my_reads.fastq --un my_reads_rRNA_unalign.fq >/dev/null

The parameters above are as follows:

   =================   ======================================================================
   Parameter           Behavior
   -----------------   ----------------------------------------------------------------------
   ``-p4``             Use 4 processors to speed alignment

   ``-v3``             Be generous, allowing up to 3 mismatches to the reeference sequence

   ``--un filename``   Save reads that don't align to `filename`
   =================   ======================================================================

For further details on alignment options, see the `bowtie`_ manual.


Chromosomal alignment
.....................

Having filtered out reads derived from rRNA, we can align the remaining data to
the genome. In this example, we use `TopHat`_, because it is capable of
performing spliced (and other gapped) alignments. We'll use splice junctions
from the ``juncs`` file made above in :ref:`starting-out-derivative-files`.

.. code-block:: shell

   # run tophat
   $ tophat -o my_reads --bowtie1 --raw-juncs splice_junctions.juncs species_bowtie_index my_reads_rRNA_unalign.fq 

   # rename output file
   $ mv my_reads/accepted_hits.bam my_reads_aligned.bam

   # index the BAM file for use in plastid, IGV, et c
   $ samtools index my_reads_aligned.bam


.. note::

   Selection of alignment parameters is a complex topic, and beyond the scope
   of this tutorial. For more details, see the `TopHat`_ manual and other forums
   dedicated to this purpose. It is important to choose alignment parameters
   tailored to your own needs.




Analysis
--------

Exploratory data analysis
.........................

TODO


Gene expression / read counting / differential expression analysis
..................................................................

See :doc:`/examples/gene_expression`



.. Metagene profiles
   .................
   TODO


Special considerations for ribosome profiling
.............................................

If analyzing :term:`ribosome profiling` data, it is helpful to go through the
following steps to assess the quality of your data.


P-site estimation
"""""""""""""""""

Estimation of ribosomal :term:`P-site offsets` is important for position-wise
analysis in ribosomal profiling. See :doc:`/examples/p_site` for background and
instructions.


Sub-codon phasing
"""""""""""""""""

After determining :term:`P-site offests`, it is possible to examine the
:term:`sub-codon phasing` found in your ribosome profiling data. See
:doc:`/examples/phasing` for information on the |phase_by_size| script.


Visualization in genome browsers
--------------------------------

Modern genome browsers can import `BAM`_ files of alignments directly. By
default, they tend to plot individual alignments and a summary track of read
coverage over each nucleotide.

Frequently it is useful to plot some aspect of the data, rather than raw read
coverage, as a function of nucleotide position in the genome. For example, in
:term:`ribosome profiling` data, this might be the total number of ribosomes
inferred to be translating a particular nucleotide. This data can be extracted
from read alignments by use of :term:`mapping rules`. See
:doc:`/concepts/mapping_rules` for a further discussion which
:term:`mapping rules` are available, and how to use them.

To export transformed data as a browser track in `bedGraph`_ or `wiggle`_
formats, see the |make_wiggle| script. If working on a large (e.g. plant or
metazoan) genome, it might be helpful to convert the `bedGraph`_ or `wiggle`_
track to a `BigWig`_ file, using `Jim Kent's utilities`_ from UCSC.


