Creating custom `BED`_, `GTF2`_, and `GFF3`_ files
==================================================

 .. TODO : update this document when custom BED columns are supported

In this tutorial, we describe how to make custom :term:`genome annotations <annotation>`
in `BED`_ and `GTF2`_ formats. These can then be used as any other annotation file:

  - to view custom regions of interest in a :term:`genome browser` like `IGV`_ or `UCSC`_
  - to perform analyses -- using pipeliness in :data:`yeti` or other programs -- on
    :term:`genomic features <feature>` that were manually constructed (e.g. reporter genes)
    or annotated

 
 .. TODO : add coordinates to glossary

Workflow
--------
While `BED`_ and `GTF2`_ are both tab-delimited formats, and could therefore be
made manually in a text editor or spreadsheet, an easier and frequently less
error-prone process is to:

 #. Enter an interactive Python session
 
 #. Manually construct |SegmentChains| and/or |Transcripts| that represent
    the :term:`features <feature>` of interest, using :term:`0-indexed`,
    :term:`half-open` coordinates
    
 #. Use the :meth:`~yeti.genomics.roitools.SegmentChain.as_bed` and
    :meth:`yeti.genomics.roitools.SegmentChain.as_gtf` methods to export
    the |SegmentChains| or |Transcripts| in the desired format

This way, all coordinate conversions, splicing, text formatting, and attribute
conversion is handled for you.


Examples
--------

 .. TODO: write examples

Reporter construct
..................
Suppose we cloned a plasmid with a 


Newly-discovered genes
......................


Choosing a format: `BED`_ vs `GTF2`_ vs `GFF3`_
-----------------------------------------------

 .. TODO: write file format choice section


Making `BigBed`_ files
----------------------
`BigBed` files are easily made from `BED`_ files using `Jim Kent's utilities`_.
To make a `BigBed`_ file:

 #. Create a custom `BED`_ file, following the examples above

 #. Sort the `BED`_ file by chromosome and start position. This is easily 
    done in a terminal session:
    
     .. code-block:: shell

        $ sort -k1,1n -k2,2n my_annotation.bed >my_annotation_sorted.bed

 #. Download and install `Jim Kent's utilities`_, which include the
    ``bedToBigBed`` program.

 #. Obtain a chromosome/contig ``.sizes`` file. If using genome builds from
    `UCSC`_, these can be downloaded using the ``fetchChromSizes`` program
    included with `Jim Kent's utilities`_.

 #. Run ``bedToBigBed``. From the terminal:

     .. code-block:: shell

        $ bedToBigBed my_annotation_sorted.bed my_genome.sizes my_annotation.bb

    Your annotation will be saved as ``my_annotation.bb``.


For more details, see the documentation for `Jim Kent's utilities`_ and the
`UCSC file format FAQ`_.

-------------------------------------------------------------------------------


See also
--------
  - :class:`~yeti.genomics.roitools.SegmentChain` and
    :class:`~yeti.genomics.roitools.Transcript` for details on these classes
  - The `UCSC file format FAQ`_ for details on file formats.
  - :doc:`/concepts/coordinates` for information on genomic coordinates

