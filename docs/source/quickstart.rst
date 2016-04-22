Getting started
===============

Genomic analysis requires some setup. This page describes the pieces that
suffice in many cases.

For those looking to try :data:`plastid` out, or to explore sequencing concepts,
we have included a :doc:`/test_dataset`, which includes sequence and annotation
for the hCMV genome, and :term:`ribosome profiling` and RNA-seq datasets. 

To get started, you need:

.. contents::
   :local:

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


Aligned sequence data
---------------------

The starting point for analysis with :data:`Plastid` is aligned sequence data,
preferably in `BAM`_ format.
 
For help on performing alignments, and a discussion of the subtleties of
choosing alignment parameters, see the documentation for the read alignment
program you use (e.g. `Bowtie`_, `Bowtie 2`_, `Tophat`_, `bwa`_). 

 
Background
----------
Most of the :data:`plastid` documentation assumes familiarty with a handful 
of concepts and conventions. We encourage those new to sequencing analysis
to check :ref:`tutorials <examples-concepts>` and browse as needed.
  

  
