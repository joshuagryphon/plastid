Creating custom BED, BigBed, and GTF2 annotation files
======================================================

In this tutorial, we describe how to make custom :term:`genome annotations <annotation>`
in `BED`_, `BigBed`_ and `GTF2`_ formats. These can then be used like any other
annotation file, for example:

 - to view custom regions of interest in a :term:`genome browser`
 
 - to perform analyses, using pipelines in :data:`plastid` or other programs

.. note::

   :data:`plastid` can also make custom `GFF3`_ files. For this see :doc:`/concepts/gff3`

   For a detailed comparison of the various annotation formats, see
   :ref:`data-annotation-format`.

This document contains the following sections:

.. contents::
   :local:

Workflow
--------
Annotation files are tab-delimited formats, and could manually be made in a text
editor or spreadsheet. However, an easier and potentially less error-prone
process is to:

#. Enter an interactive Python session

#. Construct |SegmentChains| and/or |Transcripts| that represent
   the :term:`features <feature>` of interest, using :term:`0-indexed`,
   :term:`half-open` coordinates. This can be done completely *de novo*, or 
   using other features as starting points.
   
#. Use the :meth:`~plastid.genomics.roitools.SegmentChain.as_bed` or
   :meth:`~plastid.genomics.roitools.SegmentChain.as_gtf` methods to export
   the |SegmentChains| or |Transcripts| in the desired format. To convert the
   `BED`_ file to a `BigBed`_ file, see :ref:`make-annotation-bigbed`.

:data:`plastid` will automatically handle all coordinate conversions, splicing,
text formatting, and attribute conversion.


Examples
--------

Make a GTF2 for reporter construct
..................................

Suppose we made a cell line containing a reporter construct, and that we 
want to measure the expression level of this reporter along with the expression
levels of all other genes in an :term:`RNA-seq` experiment. Let's call the vector
`lenti223` and suppose it contains a GFP and an mCherry, each driven by
some promoter.

An analysis pipeline for this experiment might look the following:

 #. Create a custom `BED`_ or `GTF2`_ file describing the reporter construct(s)
 
 #. Combine the custom annotation file with an annotation from the host genome
 
 #. Download the appropriate genome sequence, and add the vector sequence
    as another contig
    
 #. Align sequencing data to the combined genome
 
 #. Analyze sequencing data using custom annotation file.

Here, we'll focus on just making the custom `GTF2`_ file. Interactively we'll
represent the reporter transcripts as |Transcripts| and define coordinates
manually:

.. code-block:: python

   >>> from plastid import GenomicSegment, SegmentChain, Transcript

   # GFP transcript, containing 100 bp of 5' UTR and 150 bp of 3' UTR
   # 714bp coding region from bases 945-1659
   >>> gfp = Transcript(GenomicSegment("lenti223",845,1809,"+"),ID="sfGFP",cds_genome_start=945,cds_genome_end=1659)

   # mCherry transcript, similarly constructed
   >>> rfp = Transcript(GenomicSegment("lenti223",2100,3061,"+"),ID="mCherry",cds_genome_start=2200,cds_genome_end=2911)

   # now, write out features
   # we could have made a BED file using as_bed() in place of as_gtf()
   >>> with open("custom.gtf","w") as fout:
   >>>     fout.write(gfp.as_gtf())
   >>>     fout.write(rfp.as_gtf())
   >>>     fout.close()

The file ``custom.gtf`` should look something like this:

.. code-block:: shell

   lenti223    .    exon           846     1809    .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; ID "sfGFP";
   lenti223    .    CDS            946     1656    .    +    0    gene_id "gene_sfGFP"; transcript_id "sfGFP"; ID "sfGFP";
   lenti223    .    start_codon    946     948     .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; cds_start "100"; cds_end "814"; ID "sfGFP";
   lenti223    .    stop_codon     1657    1659    .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; cds_start "100"; cds_end "814"; ID "sfGFP";
   lenti223    .    exon           2101    3061    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; ID "mCherry";
   lenti223    .    CDS            2201    2908    .    +    0    gene_id "gene_mCherry"; transcript_id "mCherry"; ID "mCherry";
   lenti223    .    start_codon    2201    2203    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; cds_start "100"; cds_end "811"; ID "mCherry";
   lenti223    .    stop_codon     2909    2911    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; cds_start "100"; cds_end "811"; ID "mCherry";


Make a BED file containing halves of each coding region
.......................................................

Manually entering coordinates is laborious. More frequently, novel annotations
are derived from existing ones. Let's suppose we'd like to make a `BED`_ file
containing halves of coding regions. For this we'll use the
:doc:`demo dataset </test_dataset>`.

We'll load the transcripts, create new |SegmentChains| from those, and save
them:

.. code-block:: python

   >>> from plastid import BED_Reader

   # read transcripts   
   >>> reader = BED_Reader("some_file.bed")

   # open file for writing
   >>> halfbed = open("cds_halves.bed","w") 
   
   >>> for transcript in reader:
   >>>     cds = transcript.get_cds()
   >>>
   >>>     # make sure transcript is coding before divide CDS
   >>>     if cds.length > 0:
   >>>         name = transcript.get_name()
   >>>         halflength = cds.length // 2
   >>>
   >>>         # get halves, name each half after the parent CDS
   >>>         first_half  = cds.get_subchain(0,halflength,ID=name + "_firsthalf")
   >>>         second_half = cds.get_subchain(halflength,cds.length,ID=name + "_secondhalf")
   >>>
   >>>         # save output
   >>>         halfbed.write(first_half.as_bed())
   >>>         halfbed.write(second_half.as_bed())

   # close file
   >>> halfbed.close()



.. _make-annotation-bed-xplusy:

Making custom extended BED files
--------------------------------

:term:`Extended BED <extended BED>` and `BigBed`_ files can contain extra columns,
such as a gene ID. This can be extremely useful.

To export attributes of a |SegmentChain| or |Transcript| as extra columns
in a :term:`extended BED` format, pass a list of the attribute names (from
the dictionary `attr`) to the `extra_columns` keyword of
:meth:`SegmentChain.as_bed <plastid.genomics.roitools.SegmentChain.as_bed>`. Attributes will be
exported in the order they appear in `extra_columns`, and will be given an empty
value of "" when they are not defined

.. code-block:: python

   >>> attr = { "ID" : "some feature ID",
                "extra_field_1" : 542,
                "extra_field_2" : "some extra field",
              }

   >>> my_chain = Transcript(GenomicSegment("chrA",100,150,"+"),
                             GenomicSegment("chrA",500,550,"+"),
                             **attr)
   >>> print(my_chain.as_bed()
   chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,

   >>> print(my_chain.as_bed(extra_columns=["extra_field_1","extra_field_2"]))
   chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,    542    some extra field

If an attribute is not defined, the column will be left empty "":

.. code-block:: python

   >>> print(my_chain.as_bed(extra_columns=["extra_field_1","nonexistent_field","extra_field_2"]))
   chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,    542        some extra field



.. _make-annotation-bigbed:

Making `BigBed`_ files
----------------------
`BigBed`_ files are easily made from `BED`_ or :term:`Extended BED <extended BED>`
files using `Jim Kent's utilities`_. To make a `BigBed`_ file:

#. Create a custom `BED`_ or :term:`extended BED` file.
   For :term:`extended BED` files, consider making an optional `autoSql`_ description
   of the names & data types of the extra columns. This will allow parsers to 
   convert these to native types when reading the `BigBed`_ file.

#. Sort the `BED`_ file by chromosome and start position. This is easily 
   done in a terminal session:
   
   .. code-block:: shell

      $ sort -k1,1n -k2,2n my_annotation.bed >my_annotation_sorted.bed

#. Download and install `Jim Kent's utilities`_, which include the
   ``bedToBigBed`` program.

#. Obtain a chromosome/contig ``.sizes`` file. If using genome builds from
   `UCSC`_, these can be downloaded using the ``fetchChromSizes`` program
   included with `Jim Kent's utilities`_. For example:

   .. code-block:: shell

      $ fetchChromSizes hg38 >>hg38.sizes 

#. Run ``bedToBigBed``. From the terminal:

   .. code-block:: shell

      $ bedToBigBed my_annotation_sorted.bed my_genome.sizes my_annotation.bb

   Your annotation will be saved as ``my_annotation.bb``.

.. TODO :: add sample autoSql


-------------------------------------------------------------------------------


See also
--------
 - :ref:`data-annotation-format` for a brief overview of the costs & benefits
   of `BED`_, `BigBed`_, `GTF2`_ and `GFF3`_ files.

 - :doc:`/concepts/gff3` for information on making `GFF3`_ files
   
 - :class:`~plastid.genomics.roitools.SegmentChain` and
   :class:`~plastid.genomics.roitools.Transcript` for details on these classes
   
 - The `UCSC file format FAQ`_ for details on file formats and further discussion
   of their capabilities, advantages, and disadvantages
   
 - The `GFF3 specification <GFF3>`_ for details on GFF3 files
 
 - :doc:`/concepts/coordinates` for information on genomic coordinates
 
 - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_,
   for a description of a common `GFF3`_ feature ontology
   
 - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_,
   for the current SO consortium release.
   
 - `Jim Kent's utilities`_ for more info on making `BigBed`_ files.

 - `autoSql Grammar specification <https://github.com/ENCODE-DCC/kentUtils/blob/36d6274459f644d5400843b8fa097b380b8f7867/src/hg/autoSql/autoSql.doc>`_

