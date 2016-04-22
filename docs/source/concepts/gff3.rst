Working with GFF files
======================

GFF and its descendants are the most flexible feature-rich annotation formats.
From their flexibility arises considerable complexity. This can make them tricky
to parse. This document describes how to work with one of GFF's descendants,
`GFF3`_.


.. contents::
   :local:


Background
----------

In the original GFF specification, much of the file structure was deliberately
left unspecified to maximize flexibility. For example, it was up to users to
decide:

 - whether coordinates were :term:`0-indexed` or :term:`1-indexed`,
   and :term:`fully-closed` or :term:`half-open`
 - what feature types can be represented
 - what parent-child relationships exist (if any) between feature types
 - how the ninth column, which can contain arbitrary text data,
   including encoded attributes and values, should be formatted

Because of this, several flavors of GFF exist today, each adding
their own specifications. Two of the most common are:

 - `GTF2`_, exlusively used for describing transcripts
 - `GFF3`_, which allows a hierarchical set of feature relationships,
   and which can describe any feature ontology.

GFF3 vs GFF
-----------

`GFF3`_ adds additional constraints to the original GFF format:


 - coordinates are :term:`1-indexed` and :term:`fully-closed`
 - a ninth column contains key-value pairs of arbitary attributes
 - features can have parent-child relationships, specified by the
   `Parent` attribute in the ninth column. All features that have
   children additionally must have an `ID` attribute defined
   in the ninth column. Values of a feature's `Parent` attribute must
   match the `ID` of the parent feature.

For more detail, see the `GFF3 specification <GFF3>`_.


Parent-child relationships
--------------------------
Features in `GFF3`_ files can be hierarchical, in that they can have
*parents* and *children* defined in their ninth column. Usually, features of a
given type can only have children of specified types. For example, an `mRNA`
feature can have children of type `exon` and `CDS`, but not `gene`. A
`gene` feature, however, could have an `mRNA` child.

In addition, skipping levels of hierarchy is permitted. So,
a `gene` feature could have direct children that are `exon` or `CDS`
types, and an `mRNA` feature might not be explicitly represented.

A specification of accepted parent-child relationships is called a
*feature ontology*. The `GFF3`_ standard deliberately does not specify which
ontology to use, which adds flexibility to the format but complicates handling.


.. _gff3-feature-relationships:

Representation of discontinuous features (e.g. multi-exon transcripts)
----------------------------------------------------------------------
Only continuous features (for example exons, but not multi-exon transcripts)
can be represented directly in GFF formats.

Discontinuous features, like transcripts or gapped alignments, must be
represented as a set of related continuous components. In `GFF3`_, relationships
between a feature and its components can be represented in two ways:

 - by a shared value in the `Parent` attribute
 - by sharing an `ID` attribute

For example, a transcript can be represented as a group of exons and
coding regions, whose `Parent` attributes match the `ID` of the complex
feature::

   TODO example

Or, a transcript can be represented as a group of exons who share an `ID`::

   TODO example
    
Some research groups prefer one representation over another, while others
use both, even in the same file. 


.. _gff3-reading-overview:

Reading `GFF3`_ files in :data:`plastid`
----------------------------------------

:data:`plastid` offers two ways to read `GFF3`_ files:

 - reading them line-by-line, yielding each continuous feature as a separate
   item (via |GFF3_Reader|)
 
 - reading them processively, reconstructing transcripts from their constituent
   exons and coding regions (via |GFF3_TranscriptAssembler|)


GFF3_Reader reads individual features
.....................................
|GFF3_Reader| parses each line of a `GFF3`_ file and returns a single-segment
|SegmentChain| corresponding to the feature described by the line::
    
   >>> reader = GFF3_Reader(open("some_file.gff"))
   >>> for feature in reader:
   >>>     pass #do_something

Attributes described in the ninth column of the `GFF3`_ file are placed 
into the `attr` dictionary of the |SegmentChain|::

   >>> feature.attr
   { "some_key" : "some_value", "some_other_key" : "some_other_value", ... }


GFF3_TranscriptAssembler assembles transcripts from their component features
............................................................................

Reconstructing transcripts from `GFF3`_ files is tricky because:

 - :ref:`relationships can be represented by common Parents or shared IDs <gff3-feature-relationships>`

 - `GFF3`_ allows any possible :term:`ontology` to be used, so the relationships
   between parent and child feature types is not always clear


|GFF3_TranscriptAssembler| takes care of these two problems by:

 - first attempting to assemble transcripts by matching the `Parent` attributes
   of their component exons and coding regions, and then attempting a match
   by shared `ID` attributes
   
 - assuming that the `GFF3`_ file follows version 2.5.3 of the ontology defined
   by the `Sequence Ontology Project`_. This ontology is used by many of model
   organism databases, including `SGD`_, `FlyBase`_, and `WormBase`_.


The assembler behaves as an iterator, which assembles groups of transcripts lazily::

   >>> reader = GFF3_TranscriptAssembler("some_file.gff")
   >>> for transcript in reader: # transcripts are assembled from features when necessary
   >>>     pass # do something with the current transcript

Any malformed/unparsable `GFF3`_ lines are kept in the `rejected` property:

.. code-block:: python

   >>> reader.rejected
   [] # list of strings, corresponding to bad GFF3 lines


.. _gff3-assembly-consequences:

Memory considerations
.....................

A `GFF3`_ assembler must keep many subfeatures in memory until it is sure
that it has parsed all of the components necessary to reconstruct a given
transcript. This guarantee can be made by any of the following signals:

 - The special `GFF3`_ line `###`, which indicates that Parent-child
   relationships for preceding features have been fully resolved.
   
 - In a sorted `GFF3`_ file, a change in chromosome (assuming no feature
   spans multiple chromosomes)
   
 - The end of the `GFF3`_ file

at which point, the assembler can purge feature components from memory and
return a batch of transcripts.

Because potentially many features must be held in memory before any transcripts
can be returned, assembling transcripts from `GFF3`_ files can require
substantially more memory than reading the same data represented as 
pre-assembled transcripts line-by-line from a `BED`_ or `BigBed`_ file.

For this reason, :data:`plastid` includes a script called
:mod:`~plastid.bin.reformat_transcripts`, which can interconvert a number of
annotation formats, including `GFF3`_. 


-------------------------------------------------------------------------------

See also
--------

 - The `GFF3 specification <GFF3>`_ for a full description of the file format
 - The `Sequence Ontology Project`_ feature schema
 - |GFF3_Reader|, |GTF2_Reader|,  |GFF3_TranscriptAssembler|, and |GTF2_TranscriptAssembler|
 - :mod:`~plastid.bin.reformat_transcripts` script
