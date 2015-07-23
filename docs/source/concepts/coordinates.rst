Coordinate systems used in genomics
===================================

Genomic coordinates are typically specified as a set of:
  
  - a chromosome name
  - a start position
  - an end position
  - a chromosome strand: ('+' for the forward strand, '-' for the reverse
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
Coordinate systems can start counting from 0 (i.e. are *0-indexed*) or
from 1 (*1-indexed*). Suppose we have an XbaI restriction site on chromosome `chrI`::

                               XbaI
                              ______ 
    ChrI:         ACCGATGCTAGCTCTAGACTACATCTACTCCGTCGTCTAGCATGATGCTAGCTGAC
                  |          |^^^^^^     |          |          |
    0-index:      0          10          20         30         40 
    1-index:      1          11          21         31         41

In 0-indexed representation, the restriction site begins at position 11. In 
1-indexed representation, it begins at position 12.


In the context of genomics, both 0- and 1-indexed systems are used, depending
upon file format. :data:`yeti` knows which file formats use which representation,
and automatically converts all coordinates to a 0-indexed system, as is common
practice in `Python`_.


  .. _coordinates-half-open-fully-closed:

Half-open vs fully-closed coordinates
-------------------------------------

Similarly, coordinate systems can represent end coordinates in two ways:
 
 #. In a *fully-closed* coordinate system, positions are inclusive:
    the end coordinate corresponds to the last position **IN** the feature.
    So, in 0-indexed, fully-closed representation, the XbaI site would start at
    position 11, and end at position 16::

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

 #. In contrast, in  a *half-open* coordinate system, the end coordinate is defined as the
    first position **NOT** included in the feature. In a 0-indexed, half-open
    representation, the XbaI site woudl start at position 11, and end at 
    position 17. In this case, the length of the feature equals:

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


Following `Python`_ conventions, :data:`yeti` reports all coordinates in
0-indexed and half-open representations. In this case, the coordinate would be::

    chromosome/contig:  'ChrI'
    start:              11
    end:                17
    strand:             '.' 

    