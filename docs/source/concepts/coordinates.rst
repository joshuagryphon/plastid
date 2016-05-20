Coordinate systems used in genomics
===================================

:data:`plastid's <plastid>` readers automatically convert coordinates from 
any of the supported file formats into a :term:`0-indexed` and :term:`half-open`
space (i.e. following typical Python convention), so users don't need to worry
about off-by-one errors in their annotations.

Nonetheless, this tutorial describes various coordinate representations used
in genomics:


.. contents::
   :local:

Coordinates
-----------

Genomic coordinates are typically specified as a set of:
  
 - a chromosome name
 - a start position
 - an end position
 - a chromosome strand:
  
    - '+' for the forward strand
    - '-' for the reverse stranded
    - '.' for both strands / unstranded features

This gives rise to several non-obvious considerations:

.. _coordinates-start-end:

start ≤ end
-----------
In the vast majority of :term:`annotation` formats, the `start` coordinate
refers to the lowest-numbered (i.e. leftmost, chromosome-wise) coordinate
relative to the genome rather than the feature. So, for reverse-stand features,
the `start` coordinate actually denotes the 3' end of the feature, while the `end`
coordinate denotes the 5' end.


.. _coordinates-index-0-vs-1:

Counting from 0 vs 1
--------------------
Coordinate systems can start counting from 0 (i.e. are :term:`0-indexed`) or
from 1 (:term:`1-indexed`). Suppose we have an XbaI restriction site on
chromosome `chrI`:

.. code-block:: none

                              XbaI
                             ______ 
   ChrI:         ACCGATGCTAGCTCTAGACTACATCTACTCCGTCGTCTAGCATGATGCTAGCTGAC
                 |          |^^^^^^     |          |          |
   0-index:      0          10          20         30         40 
   1-index:      1          11          21         31         41

  

In :term:`0-indexed` representation, the restriction site begins at position 11.
In :term:`1-indexed` representation, it begins at position 12.

In the context of genomics, both :term:`0-indexed` and :term:`1-indexed`
systems are used, depending upon file format. :data:`plastid` knows which file
formats use which representation, and automatically converts all coordinates
to a :term:`0-indexed` representation, following Python idioms.


.. _coordinates-half-open-fully-closed:

Half-open vs fully-closed coordinates
-------------------------------------

Similarly, coordinate systems can represent end coordinates in two ways:
 
#. In a :term:`fully-closed` or :term:`end-inclusive` coordinate system,
   positions are inclusive: the end coordinate corresponds to the last
   position **IN** the feature.

   So, in :term:`0-indexed`, :term:`fully-closed` representation,
   the XbaI site would start at position 11, and end at position 16::

                                XbaI
                               ______ 
     ChrI:         ACCGATGCTAGCTCTAGACTACATCTACTCCGTCGTCTAGCATGATGCTAGCTGAC
                   |           ^^^^^^     |          |          |
     0-index:      0           |    |     20         30         40 
                               |    |
     Start & end:              11   16
                                 
   And the length of the feature equals:

   .. math::
     
       \ell = end - start + 1 = 16 - 11 + 1 = 6

#. In contrast, in  a :term:`half-open` coordinate system, the end coordinate
   is defined as the
   first position **NOT** included in the feature. In a :term:`0-indexed`,
   :term:`half-open` representation, the XbaI site starts at position 11, and
   ends at position 17. In this case, the length of the feature equals:

   .. math::
     
       \ell = end - start = 17 - 11 = 6


Four possible coordinate representations
----------------------------------------
Because coordinate systems can be :ref:`0-indexed or 1-indexed <coordinates-index-0-vs-1>`,
and :ref:`half-open or fully-closed <coordinates-half-open-fully-closed>`,
genomic features can be can be represented in four possible ways. For the XbaI
site in this example:

   =============   =============    ==================
        \          **Half-open**    **Fully-closed**
   -------------   -------------    ------------------
   **0-indexed**   start: 11        start: 11
                   end: 17          end: 16

   **1-indexed**   start: 12        start: 12
                   end: 18          end: 17

   =============   =============    ==================


Coordinate systems of some common file formats
----------------------------------------------

   =============   =============   ====================
   **Format**      **Index**       **End coordinates**
   -------------   -------------   --------------------
   `BED`_          0               Half-open
   `BigBed`_       0               Half-open
   `GTF2`_         1               Fully-closed
   `GFF3`_         1               Fully closed
   Other GFFs      Either          Either
   `PSL`_          0               Half-open
   -------------   -------------   --------------------
   `SAM <BAM>`_    1               n/a
   `BAM`_          0               n/a
   bowtie          0               n/a
   -------------   -------------   --------------------
   `bedGraph`_     0               Half-open
   `BigWig`_\*     0 or 1          Half-open or n/a          
   `Wiggle`_       1               n/a
   =============   =============   ====================
 
\*The coordinate representation used in `BigWig`_ files depends upon
the format of the data blocks inside the file, which can be represented
as `wiggle`_ or `bedGraph`_ blocks.


Conventions used in `plastid`
-----------------------------
Following `Python`_ conventions, :data:`plastid` reports all coordinates in
:term:`0-indexed` and :term:`half-open` representation.
In this case, the coordinate would be:

.. code-block:: none

   chromosome/contig:  'ChrI'
   start:              11
   end:                17
   strand:             '.' 


-------------------------------------------------------------------------------

See also
--------
 - `UCSC file format FAQ`_ for detailed descriptions of various file formats
