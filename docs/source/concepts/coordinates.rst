Coordinate systems in genomics
==============================

Genomic coordinates are typically specified as a set of:
  
  - a chromosome name
  - a start position
  - an end position
  - a chromosome strand ('+' for the forward strand, '-' for the reverse
    strand, or '.' for a feature on both/no strands).

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
from 1 (*1-indexed*). In the context of genomics, both systems are used,
depending upon file format.

:data:`yeti` knows which formats use which system, and therefore converts
all coordinates to a 0-indexed system, as is common practice in `Python`_.


  .. _coordinates-half-open-fully-closed:

Half-open vs fully-closed coordinates
-------------------------------------
Let us suppose we have a feature that is nine nucleotides long, and begins
on chromosome I, at base 12. This means the feature includes bases
12,13,14, 15,16,17, 18,19, and 20.

We can describe its coordinates in two ways:

 #. In a *fully-closed* coordinate system, positions are inclusive. So,
    the end coordinate corresponds to the last position **IN** the feature.
    Therefore:
      - `start = 12`
      - `end = 20`
    
    And, the length of the feature equals:
    
        `end - start + 1 = 20 - 12 + 1 = 9`
 
 #. In a *half-open* coordinate system, `end` is defined as the first position
    **NOT** included in the feature. Therefore, we would write:
    
      - `start = 12`
      - `end = 21`
    
    And the length of the feature is simply:
    
        `end - start = 21 - 12 = 9`

Some :term:`annotation` formats use fully-closed coordinates, others
half-open. Again, :data:`yeti` takes care of this for you, and converts
all coordinates to half-open space, in keeping with `Python`_ conventions.


.. TODO: make graphic

Example: a restriction site
---------------------------

Suppose we have a very short genome with the following sequence. It contains
an XbaI (*TCTAGA*) restriction site, which have highlighted below::

                               XbaI
                              ______ 
    Sequence:     ACCGATGCTAGCTCTAGACTACATCTACTCCGTCGTCTAGCATGATGCTAGCTGAC
                              ^^^^^^

As mentioned :ref:`above <coordinates-index-0-vs-1>`, coordinates for this
sequence can be *0-indexed*, or *1-indexed*, which means that our restriction
site starts either at base 11 or 12, respectively::


                               XbaI
                              ______ 
    Sequence:     ACCGATGCTAGCTCTAGACTACATCTACTCCGTCGTCTAGCATGATGCTAGCTGAC
                  |          |^^^^^^     |          |          |
    0-index:      0          10          20         30         40 
    1-index:      1          11          21         31         41

Coordinate systems can also be :ref:`half-open or fully-closed <coordinates-half-open-fully-closed>`.
Therefore, the restriction site can be described in four possible ways:

    =============   =============    ==================
         \          **Half-open**    **Fully-closed**
    -------------   -------------    ------------------
    **0-indexed**   start: 11        start: 11
                    end: 17          end: 16

    **1-indexed**   start: 12        start: 12
                    end: 18          end: 17
    =============   =============    ==================


:data:`yeti` is 0-indexed and half-open, and would therefore report the coordinates
as::

    chromosome/contig:  'sequence'
    start:              11
    end:                17
    strand:             '.' 
