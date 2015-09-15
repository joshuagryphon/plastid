Working with GFF files
======================


GFF files are the most feature-rich annotation format. In the original
specification, much of the file strucure was deliberately left
unspecified, to allow users flexibility to change what they needed.
For example, it was up to users to decide:

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

In this tutorial, we discuss briefly discuss handling `GFF3`_ files
in :data:`plastid`. We recommend that users read the `GFF3 specification <GFF3>`_
for full discussion of the file type. In `GFF3`_:

  - coordinates are :term:`1-indexed` and :term:`fully-closed`
  - a ninth column containing key-value pairs of arbitary attributes
  - that features can have parent-child relationships, specified by the
    `Parent` attribute in the ninth column. All features that have
    children additionally must have an `ID` attribute defined
    in the ninth column. Values of a feature's `Parent` attribute must
    match the `ID` of the parent feature.


Parent-child relationships
--------------------------
Features in `GFF3`_ files can be hierarchical, in that they can have
parents and children. Usually, features of a given type can only
have children of specified types. For example, an `mRNA` feature
can have children of type `exon` and `CDS`, but not `gene`. A
`gene` feature, however, could have an `mRNA` child.

In addition, skipping levels of hierarchy is permitted. So,
a `gene` feature could have direct children that are `exon` or `CDS`
types, and an `mRNA` feature might not be explicitly represented.

A specification of accepted parent-child relationships is called
a feature :term:`ontology`. The `GFF3`_ standard deliberately does
not specify which :term:`ontology` to use, which adds flexibility
to the format but can complicate analysis.


 .. _gff3-feature-relationships

Complex features
----------------
Discontinuous features, like transcripts or gapped alignments, can be
represented in `GFF3`_ files as a set of related continuous sub-features.
Relationships between sub-features can be represented in two ways:

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
-------------------------------------

Reading simple features
.......................
|GFF3_Reader| parses each line of a `GFF3`_ file and returns a single-segment
|SegmentChain| corresponding to the feature described by the line::
    
    >>> TODO open GFF3

Attributes described in the ninth column of the `GFF3`_ file are placed 
into the `attr` dictionary of the |SegmentChain|::

    >>> feature.attr


Assembling complex features
...........................
One reason to use a `GFF3`_ file is to preserve relationships between features,
and/or to assemble complex, discontinuous features.

However, because :ref:`relationships can be represented by common Parents or shared IDs <gff3-feature-relationships>`,
and because `GFF3`_ is agnostic to the feature  :term:`ontology` used, 
correctly assembling complex features from a `GFF3`_ file is not trivial.

For convenience, |GFF3_TranscriptAssembler| is provided.
By assuming the `GFF3`_ uses th SO 2.5.3 ontology (used by many of the model
organism databases,  including `SGD`_, `FlyBase`_, and `WormBase`_), it 
assembles features into |Transcript| objects, first by `Parent` matching,
and then by shared `ID`, if shared `ID` attributes are present.

The reader behaves as an iterator, which assembles groups of transcripts lazily::

    >>> reader = GFF3_TranscriptAssembler("some_file.gff")
    >>> for transcript in reader:
    >>>     pass # do something

Any malformed/unparsable `GFF3`_ lines are kept in the `rejected`
attribute::

    >>> reader.rejected
    [] # list of strings, corresponding to bad GFF3 lines


 .. _gff3-assembly-consequences:

Consequences of assembly
........................
Because complex features are made of sub-features, a `GFF3`_ assembler 
must keep many features in memory until it is certain it has collected
all of a feature's sub-features. Three signals can indicate when
it is time to assemble:

  - The special line `###`, which indicates that Parent-child relationships
    for preceding features have been fully resolved.
  - In a sorted `GFF3`_ file, a change in chromosome (assuming no feature
    spans multiple chromosomes)
  - The end of the `GFF3`_ file

Because so many features must be held in memory before a feature can 
be assembled from subfeatures, assembling a transcript form a `GFF3`_ file
requires much more memory than simply reading a transcript from a single
line of a `BED`_ file.


 .. _gff3-write-assembler: 

Writing your own assembler
..........................
It is possible to write custom assemblers transcripts (or any complex feature)
from any :term:`ontology`. |AbstractGFF_Assembler| is provided 
as a base class.

handling attributes? pooled attribute func
stop feature?

TODO : finish section on writing own assembler



-------------------------------------------------------------------------------

See also
--------

  - The `GFF3 specification <GFF3>`_ for a full description of the file format
  - The Sequence Ontology consortium's feature schema
  - |GFF3_Reader|, |GTF2_Reader|,  |GFF3_TranscriptAssembler|, and |GTF2_TranscriptAssembler|

