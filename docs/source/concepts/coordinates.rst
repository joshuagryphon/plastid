Coordinate systems in genomics
==============================
For genomic entities, coordinatesare typically specified as a set of:
  
  - a chromosome name
  - a start position
  - an end position
  - a chromosome strand ('+' for the forward strand, '-' for the reverse
    strand, or '.' for a feature on both/no strands).

This seems clear enough, but there are several non-obvious considerations:

`start <= end`
--------------
In the vast majority of :term:`annotation` formats, the `start` coordinate
refers to the lowest-numbered (i.e. leftmost, chromosome-wise) coordinate,
rather than the fiveprime end of a feature. So, for reverse-stand features,
the `start` coordinate actually denotes the feature's 3' end.

Counting from 0 vs 1
--------------------
Coordinate systems typicaly start counting from 0 (i.e. are *0-indexed*) or
from 1 (*1-indexed*). In the context of genomics, both systems are used,
depending upon file format. In this case, some formats label the first
nucleotide of a chromosome `0`, while others label the first nucleotide `1`.

:data:`yeti` knows which formats use which system, and therefore converts
all coordinates to a 0-indexed system, as is common practice in `Python`_.

Half-open vs fully-closed coordinates
-------------------------------------
Let us suppose we have a feature that is nine nucleotides long, and begins
on chromosome I, at base 12. This means the feature in includes bases
12,13,14, 15,16,17, 18,19, and 20.

We can describe its coordinates in two ways:

 #. In a fully-closed coordinate system, positions are inclusive. So,
    the end coordinate corresponds to the last position *IN* the feature.
    Therefore:
      - `start = 12`
      - `end = 20`
    
    And, the length of the feature equals:
    
        `end - start + 1 = 20 - 12 + 1 = 9`
 
 #. In a half-open coordinate system, `end` is defined as the first position
    *NOT* included in the feature. Therefore, we would write:
    
      - `start = 12`
      - `end = 21`
    
    And the length of the feature is simply:
    
        `end - start = 21 - 12 = 9`

Some :term:`annotation` formats use fully-closed coordinates, others
half-open. Again, :data:`yeti` takes care of this for you, and converts
all coordinates to half-open space, in keeping with `Python`_ conventions.

