Getting started
===============

Genomic analysis requires some setup. This page provides a quick overview
of those pieces.


To get started, you need:

.. contents::
   :local:


For those looking to try :data:`plastid` out, or to explore sequencing concepts,
we have included a :doc:`/test_dataset`, which includes sequence and annotation
for the hCMV genome, and :term:`ribosome profiling` and RNA-seq datasets. 

For those setting up their own data, please continue reading:


A genome sequence & annotation
------------------------------

The starting point for most genomics research is to obtain a genome sequence
and matching :term:`genome annotation <annotation>`. Good sources for these
include:

  - Mammal genomes: `UCSC`_, `Ensembl`_, `ENCODE`_, and `GENCODE`_
  - Fly genomes: `FlyBase`_, `modENCODE`_
  - *S. cerevisiae*: `SGD`_

It is critical that the genome sequence and feature annotations use the same
coordinates, so be sure to download corresponding versions from a single build
(i.e. it is unhelpful to mix mouse the *mm9* genome sequence with the *mm10*
annotation).

Often it is useful to do some pre-processing of files once they have been
downloaded. Detailed discussion is provided in :doc:`/examples/a1_genome_setup`


Aligned sequence data
---------------------

The starting point for analysis with :data:`Plastid` is aligned sequence data,
preferably in `BAM`_ format.

An brief overview of the relevant steps in setting up alignments and exploring 
data may be found in  :doc:`/examples/a2_sample_processing`.

That said, choice of alignment parameters merits careful consideration, which
is a weighty topic, beyond the scope of this tutorial. For a more detailed
discussion, see the documentation for the read alignment program you use (e.g.
`Bowtie`_, `Bowtie 2`_, `Tophat`_, `bwa`_, `STAR`_).


 
Other background info
---------------------

Most of the :data:`plastid` documentation assumes familiarty with a handful 
of concepts and conventions. We encourage those new to sequencing analysis
to check the :doc:`/examples` and :doc:`/glossary` as needed.
  

  
