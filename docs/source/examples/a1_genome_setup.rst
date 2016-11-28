Setting up a genome for analysis
================================

When on-boarding a new genome, it helps to do some pre-processing. Most of this
involves downloading and converting files into formats that can be read
more quickly. Some of it involves creating new files.

Details on why this is important can be found in :doc:`/concepts/data`.

Fortunately, these steps only need to be performed when the genome annotation or
sequence changes.

.. note::

   This tutorial is under construction. 

.. note::

   The examples in this tutorial are **sequential**. You don't need to do every 
   step, but some steps require subsequent steps to be modified if they are 
   skipped. Instructions for doing so are included in each step. 
   
   If you miss this information, you may get spurious results or frustrating errors.


This document includes the following sections:

.. contents::
   :local:


Download some tools
-------------------

In this tutorial, we use the following tools. They are useful to have on-hand
in many circumstances:

 - `SAMtools`_ and `htslib`_ to index and manipulate `BAM`_ files

 - `bowtie`_ and `Tophat`_, or your favorite short read aligners (`bwa`_,
   `star`_, `kallisto`_, et c)

 - The Hannon Lab's `fastx toolkit`_ for removing cloning adaptors
 
 - `Jim Kent's utilities`_ for making `BigBed`_ or `BigWig`_ files



.. _starting-out-annotation:
 
Obtain a genome sequence & matching annotation
----------------------------------------------

Choose a curator
................

For many organisms, multiple genome assemblies exist. Often, there are multiple
versions, and also multiple curators. For example, `Ensembl`_, `UCSC`_, and
`RefSeq`_ all host versions of the mouse and human genomes.

Because chromosomal sequences differ between assemblies, it is essential that
the genome annotation (containing the locations of genes and transcripts) come
from the same curator as assembly (sequences). It is also important that assembly
and annotation verions match.

In other words, if using the current `UCSC`_ mouse genome assembly, be sure to
use their annotation as well. Do not mix and match!

Sources for genome assemblies and annotations include:

 - `SGD`_ for *S. cerevisiae* S228C (reference strain), as well as many others
   (available in their `downloads` page)

 - `FlyBase`_ and `modENCODE`_ for *D. melanogaster* and other insects 
 
 - `ENCODE`_/`GENCODE`_, `Ensembl`_, `UCSC`_, and `RefSeq`_ for human, mouse,
   and other model organsisms
   
 - `NCBI Assembly`_ and `JGI`_ for others (e.g. non-model organisms
   alternate strains, data hot off the press)


Download files
..............

Once you have chosen a data source, you will need:

 - **Genome sequence**: this will generally be in `fasta`_ or `twobit`_ format.
   Occasionally, it will be included in a `GFF3`_ file, in which case it
   will need to be extracted.
   
 - **Genome annotation**: this will appear in `GTF2`_ (just transcripts and
   coding regions), `GFF3`_ (all types of features- e.g. genes, origens of
   replication, microRNAs, et c), or `BED`_ format.
   
   `BED`_ format files are easy to work with, but usually contain transcripts
   divorced from their genes. We recommend grabbing the `GTF2`_ or `GFF3`_ file
   as well as a `BED`_ file if available. In many cases, we'll re-make a new
   `BED`_ file using the info from the `GTF2`_/`GFF3`_ later.


.. _starting-out-derivative-files:

Convert file formats
--------------------


Extract genome sequence from GFF3
.................................

If the annotation and/or sequence is in `GFF3`_ format, we made need to
extract the genome sequence if a separate `fasta`_ file is not available.
To make a `fasta`_ file from data embedded in a `GFF3`_, open a terminal
and type:

.. code-block:: shell

   $ cat my_file.gff | awk 'BEGIN { doprint = 0}; \
                            doprint == 1 { print $0 }; \
                            $0 ~ /#FASTA/ { doprint = 1 }' >my_file.fa


Make an indexed sequence file
.............................

If you're using a large genome and plan to do lots of sequence manipulation
downstream, consider making either a `twobit`_ file or an indexed `bgzf`_-compressed
`fasta`_ file from the `fasta`_ file. `BioPython`_ reads `bgzf`_-compressed
`fasta`_ files  natively, and `twobitreader`_ reads `twobit`_ files.
:data:`plastid` works with objects produced by both packages.

To compress the `fasta`_ file with `bgzf`_, make sure to have `htslib`_ installed,
then type from the terminal:

.. code-block:: shell

   $ bgzip my_file.fa

To make a `twobit`_ file, make sure to have `Jim Kent's utilities`_ installed,
then type from the terminal:

.. code-block:: shell

   $ faToTwoBit my_file.fa my_file.2bit

..

Sort your files
...............

Because of the way features are stored in `GFF3`_  and `GTF2`_, we can get away
with using far less memory if they are sorted. In general, it is good practice
to sort annotation files by chromosome, starting coordinate, ending coordinate,
and strand.

There are multiple ways to sort these file. One easy one is to use the ``sort``
utility. Follow the examples below:

.. code-block:: shell

   # sort a GFF3 file
   $ cat my_file.gff | grep -v "#" | sort -k1,1 -k4,4n >my_file_sorted.gff
   
   # sort a GTF2 file - the same way!
   $ cat my_file.gtf | grep -v "#" | sort -k1,1 -k4,4n >my_file_sorted.gtf
   
   # sort a BED file
   $ cat my_file.bed | grep -v "#" | -k1,1 -k2,2n >my_file_sorted.bed

When working with sorted files in :data:`plastid`, use the ``--sorted`` 
command-line argument for command line scripts (as in the examples below), or
the keyword argument ``is_sorted=True`` when using various readers interactively
(for an example, see the documentation for
:class:`~plastid.readers.gff.GTF2_TranscriptAssembler`).

.. _starting-out-convert-gff3:

Assemble transcripts from GFF3 or GTF2
......................................

`GFF3`_ files include many genomic features aside from transcripts. It often
helps to filter these out, and make a more compact, more standardized `GTF2`_
or `BED`_ file.

That said, both `GFF3`_ and `GTF2`_ files are memory hogs, and difficult to use
with large (e.g. mammal, ciliate, plant) genomes. `BED`_ files use far less
memory, but in their native format fail to include useful information like
parent-child relationships between genes and transcripts.

A useful compromise is to use an :term:`extended BED` file or a `BigBed`_ file,
which are still memory efficient but can preserve extra information like
gene-to-transcript relationships. However, comparatively few tools support these 
formats (:data:`plastid` does), while `GFF3`_, `GTF2`_ and `BED`_ are more or
less universally  supported.

So, here we provide information on how to make all of them. :data:`plastid`
includes a script called |reformat_transcripts|  for this. If starting with
a `GTF2`_ file, replace ``--annotation_format GFF3`` with ``--annotation_format GTF2``
in the examples below. We assume your files are sorted by chromosome. If not,
drop the ``--sorted`` argument from below, and go do laundry or something
while things run.

To make a `GTF2`_ file from a `GFF3`_ file:

.. code-block:: shell

   $ reformat_transcripts --annotation_files my_file.gff \
                          --annotation_format GFF3 \
                          --sorted \
                          --output_format GTF2 \
                          myfile.gtf2


To make a `BED`_ file:

.. code-block:: shell

   $ reformat_transcripts --annotation_files my_file.gff \
                          --annotation_format GFF3 \
                          --sorted \
                          --output_format BED \
                          myfile.bed

To make an :term:`extended BED` file that include extra columns for
gene-transcript relationships and notes (assuming an attribute called `Notes`
is defined in the `GFF3`_/`GTF2`_ file), add the ``--extra_columns`` argument:

.. code-block:: shell

   $ reformat_transcripts --annotation_files my_file.gff \
                          --annotation_format GFF3 \
                          --output_format BED \
                          --sorted \
                          --extra_columns gene_id Notes\
                          -- \
                          myfile_extended.bed


.. _starting-out-make-bigbed:

To convert a `BED`_ or :term:`extended BED` file to a `BigBed` file, use the
``bedToBigBed`` program from `Jim Kent's utilities`_.  In this example,
we'll use the :term:`extended BED` file.

To take advantage of the extra columns we output, ``bedToBigBed`` needs to know
what they contain. Fortunately, the screen output from |reformat_transcripts|
will contain a table declaration for this purpose. In this example, it reads:

.. code-block:: none

   table bigbed_columns "myfile_extended columns"
   (
       string            chrom;          "Chromosome"
       uint              chromStart;     "chr start"
       uint              chromEnd;       "chr end"
       string            name;           "item name"
       uint              score;          "score"
       char[1]           strand;         "strand"
       uint              thickStart;     "thickstart"
       uint              thickEnd;       "thickend"
       uint              reserved;       "normally itemRgb"
       int               blockCount;     "block count"
       int[blockCount]   blockSizes;     "block sizes"
       int[blockCount]   chromStarts;    "block starts"
       string            gene_id;        "description of custom field contents"
       string            Notes;          "description of custom field contents"
   )
   
   
Copy and paste that into a file called ``my_fields.as``, editing the descriptions
of `gene_id` and `Notes` as you see fit.

Next, we need to know the chromosome sizes so the `BigBed`_ file can be 
indexed. ``bedToBigBed`` expects a two-column tab-delimited text file where
the columns are `(chromosome name, chromosome size)`. If you have the genome
sequence as a `fasta`_ file, create a file of chromosome sizes by entering
the following in a Python terminal:

.. code-block:: python

   >>> from Bio import SeqIO
   >>> genome = SeqIO.parse(open("my_genome_sequence.fa"),"fasta")
   
   >>> outfile = open("mychroms.sizes","w")
   >>> for my_seq in genome:
   >>>     outfile.write(my_seq.id + "\t" + str(len(my_seq)) + "\n")
   >>>
   
   >>> outfile.close()
   

Finally, sort the :term:`extended BED` file and run ``bedToBigBed``. If you're
using the non-extended `BED`_ file, drop the arguments
``-type=12+2 --extraIndex=name,gene_id, -as=my_fields.as``
from the example. Type from the bash terminal:

.. code-block:: shell

   $ sort -k1,1 -k2,2n myfile_extended.bed >myfile_sorted.bed
   $ bedToBigBed -tab -type=bed12+2 \
                 -extraIndex=name,gene_id \
                 -as=my_fields.as \
                 myfile_sorted.bed mychroms.sizes myfile.bb 

The output file, ``myfile.bb``, can be used in :data:`plastid` and by other
programs and genome browsers that support `BigBed`_ files. The use of
``-extraIndex`` (if using the :term:`extended BED` file) is a bonus- it throws
extra indices on the names and gene IDs
of transcripts in the `BigBed`_ file, enabling :data:`Plastid` to search the
file efficiently by name. See :mod:`plastid.readers.bigbed` for examples.


Using any of these files with :data:`plastid`'s command-line scripts is easy.
Just be sure to set the ``--annotation_format`` argument to ``BED``, ``BigBed``,
``GTF2`` or ``GFF3``, depending upon which format you are using.

If using an :term:`extended BED` file, set ``--annotation_format`` to ``BED``,
and use the ``--bed_extra_columns`` argument to specify the column names. For example,
to convert the :term:`extended BED` back to a `GTF2`_:

.. code-block:: shell

   $ reformat_transcripts --annotation_files myfile_extended.bed \
                          --annotation_format BED \
                          --bed_extra_columns gene_id Notes\
                          --sorted \
                          --output_format GTF2 \
                          myfile.gtf
   
..

Consider excluding transcripts supported by scant evidence
..........................................................

Curators of genome annotations require different thresholds of evidence when
adding transcripts to the annotation.

Some curators are very conservative, and at the risk of exluding valid
transcripts, require lots of verification to add a transcript model. Others are
very inclusive, and include many potentially dubious trancripts in their
annotations. Obviously, this is a tradeoff, and depending upon your needs,
you may wish to be conservative or inclusive.

A discussion of whether or when to filter transcripts is beyond the scope of
this tutorial. However, removing transcripts that have little or no biological
evidence can improve processing times, in particular for vertebrate genomes.

A useful resource for this is the `APPRIS`_ database (:cite:`Rodriguez2013`),
which uses multiple sources of evidence (e.g. conservation among vertebrates,
presence of protein-coding domains, et c) to label transcript isoforms
as `principal`, `alternative`, or `minor`, and within each of these categories,
provides a reliability score.

One option, then, is to use `APPRIS`_ annotations to remove all isoforms
considered `minor`, retaining those considered `principal` or `alternative`.


Then, re-sort, index, convert the resulting annotation file as described 
in :ref:`the previous section <starting-out-convert-gff3>`.



Make files for aligners
-----------------------


bowtie indexes
..............

`bowtie`_ and `Tophat`_ require pre-built genome indexes before alignment.
These are built from `fasta`_ files of genome sequence. From the terminal:

.. code-block:: shell

   $ bowtie-build my_file.fa my_genome_index
   
This will create six files, all beginning with `my_genome_index`.


.. _starting-out-juncs-file: 

.juncs file for Tophat
......................

`Tophat`_ uses a custom file format (``.juncs``) to specify splice junctions.
While `Tophat`_ can extract junctions from a `GTF2`_ file, it is often convenient
to have a pre-built ``.juncs`` file. Plastid includes a script called |findjuncs|
for this. To use it, type from the terminal:

.. code-block:: shell 

   $ findjuncs --annotation_files my_file.gtf \
               --annotation_format GTF2 \
               --sorted \
               --export_tophat \
               my_junctions
   
This will create two files: ``my_junctions.juncs`` and ``my_junctions.bed``, a
`BED`_ format file containing the same splice junctions represented as exon pairs.
|findjuncs| can perform a number of other useful pieces of information- see
its documentation to see what it can do.



Make files for analysis in Plastid
----------------------------------

If you plan to do your analysis using :data:`plastid`'s tools, it is helpful
to pre-compute a number of files, in the following order:


Crossmap
........

Repetitive elements and recently-duplicated genes often contain sequence that
can :term:`multimap`, or align equally well to multiple places in the genome.
In a sequencing dataset, the origin of multimapping reads is often ambiguous,
and it can be useful to exclude such reads, as well as the regions from the 
genome that give rise to them, from analysis. Plastid includes a tool called
|crossmap| for this purpose. The algorithm is described in detail
:ref:`here <masking-crossmap-script>`.
 
To use |crossmap|, you need a sequence file, a `bowtie`_ index, and `bowtie`_ 
(not `bowtie 2`_) installed on your machine. Follow the example below, but 
set the parameter ``-k`` to the expected size of your read alignments (or a 
little shorter, to be conservative):

.. code-block:: shell

   $ crossmap -k 26 --mismatches 2 my_file.fa /path/to/bowtie/indexes/my_genome_index my_crossmap 


For further discusison of |crossmap| and examples of how to use its output,
see :doc:`/examples/using_masks`. 

As with all annotation files (and, a :term:`crossmap` is an annotation), we
recommend sorting and, for large genomes, converting the output from `BED`_
to `BigBed`_ format. Instructions for how to do this can be found in the screen
output from |crossmap|, or  :ref:`above <starting-out-make-bigbed>`.

Examples below assume you have made a :term:`crossmap`. If you decide not to,
drop all of the ``--mask_annotation*`` arguments from the examples below.

.. _starting-out-maximal-spanning-windows: 

Maximal spanning windows
........................

:data:`plastid` uses :term:`maximal spanning windows <maximal spanning window>`
for :term:`metagene analysis <metagene>`, estimation of
:term:`P-site offsets <P-site offset>`, and determination of
:term:`sub-codon phasing` for :term:`ribosome profiling` data.

If you plan to do any of these analyses, complete these steps, use
:data:`plastid`'s |metagene| script to build
:term:`maximal spanning windows <maximal spanning window>` from a genome annotation.


The algorithm for window generation is described in depth in the
:ref:`metagene tutorial <metagene-generate>`. To generate the windows, type from
the terminal:

.. code-block:: shell

   $ metagene generate --annotation_files my_file.gtf --sorted \
              --mask_annotation_files my_crossmap.bb \
              --mask_annotation_format BigBed \
              --downstream 200 \
              my_windows
              
..


Corrected gene boundaries
.........................

Plastid includes tools for measuring gene expression. One of these, called |cs|,
includes a ``generate`` mode, that creates a corrected genome annotation. This
correct annotation differs from an off-the-shelf genome annotation in that it:

 - collapses genes whose transcripts share exact exons to "merged" genes 

 - excludes regions of the genome that are overlapped by more than one merged
   gene
 
 - excludes regions flagged in a :term:`mask file` or :term:`crossmap`
 
 - classifies sub-regions of genes as *CDS*, *5' UTR*, or *3' UTR* based upon
   how they appear in that gene's transcripts

The corrected annotation may then be used by `cs` or other tookits (within 
:data:`plastid` or outside it), for gene expression measurement. For details,
see the documentation for |cs|.




Sessions for genome browsing
----------------------------

We can't stress enough the value of visual examination of your data in a 
genome browser. We like The Broad Institute's `IGV`_ for this purpose. To
create a session for your genome, follow the instructions for
`creating a .genome file, here <http://software.broadinstitute.org/software/igv/NewGenomeMgmt>`_.

Supply your custom `GTF2`_, `BED`_, or `BigBed`_ file as the `gene file`,
and the exact genome sequence (`fasta`_ format) used above as the `sequence file`.

We also recommend loading a `BigBed`_ version of any crossmap file; it will help
explain curious increases or decreases in read density that you might see when
you load your alignments.



