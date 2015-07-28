Creating custom `BED`_, `GTF2`_, and `GFF3`_ files
==================================================

 .. TODO : update this document when custom BED columns are supported

In this tutorial, we describe how to make custom :term:`genome annotations <annotation>`
in `BED`_ and `GTF2`_ formats. These can then be used as any other annotation file:

  - to view custom regions of interest in a :term:`genome browser` like `IGV`_ or `UCSC`_
  - to perform analyses -- using pipelines in :data:`yeti` or other programs -- on
    :term:`genomic features <feature>` that were manually constructed (e.g. reporter genes)
    or annotated


Workflow
--------
Annotation files are tab-delimited formats, and could therefore be
made manually in a text editor or spreadsheet, an easier and frequently less
error-prone process is to:

 #. Enter an interactive Python session
 
 #. Manually construct |SegmentChains| and/or |Transcripts| that represent
    the :term:`features <feature>` of interest, using :term:`0-indexed`,
    :term:`half-open` coordinates
    
 #. Use the :meth:`~yeti.genomics.roitools.SegmentChain.as_bed`,
    :meth:`~yeti.genomics.roitools.SegmentChain.as_gtf`, and 
    :meth:`~yeti.genomics.roitools.SegmentChain.as_gff3` methods to export
    the |SegmentChains| or |Transcripts| in the desired format

All coordinate conversions, splicing, text formatting, and attribute
conversion is handled for you.


Examples
--------

 .. TODO: write examples

Reporter construct
..................
[TODO]


A complex feature that is not a transcript
..........................................
[TODO]
 .. Use this to explain problems with GFF3 export ontologies
 

Choosing a format: `BED`_ vs `GTF2`_ vs `GFF3`_
-----------------------------------------------

Which format to choose depends on your purposes. The following questions
highlight the differences between formats:

Are your features something other than exons, coding regions, UTRs or transcripts?
..................................................................................
If yes, you need to use something other than `GTF2`_, which, according to its
specification, can only transcripts and their parts:

  - CDS
  - start_codon
  - stop_codon
  - 5UTR
  - 3UTR
  - inter
  - inter_CNS
  - intron_CNS
  - exon
  

Do you need rich attribute data?
................................
If yes, use something other than `BED`_ format. In classic `BED`_ formats,
it is possible to store the only the following attributes:

  - feature name
  - feature coordinates (feature can be discontinuous, like a multi-exon transcript)
  - feature coding region start & stop  
  - a score for the feature
  - a color for rendering the feature in a genome browser

Any other data (e.g. GO terms, IDs of parental or related features, et c) cannot
be represented.

`BigBed`_ and `BED+X`_ formats can include additional attributes
as columns, but in these formats all records must contain the same types of
attributes.

`GTF2`_ and `GFF3`_ offer the richest feature descriptions because they contain
a specific column (column 9) that holds key-value pairs describing arbitrary
information, which can differ from record to record.


Do your dataset have multiple types of features?
................................................
If so, use `GTF2`_ or `GFF3`_. Because `BED`_ files contain no column to
describe feature type, it is simplest to make sure all features in the `BED`_
file are of a single type.


Are the features you care about discontinuous? And is your computer limited for memory?
.......................................................................................
If yes, use one of the `BED`_-family formats. In `BED`_ files, each feature
-- even discontinuous
features like multi-exon transcripts -- are represented as single lines. This
means that programs don't need to search through a file to find all of the
pieces (e.g. exons, piecse of coding regions, et c) that make up a feature.

In contrast, in `GTF2`_ and `GFF3`_ files, each line can only contain a continuous
feature or sub-feature. So, to represent a multi-exon transcript, each exon
would be represented on its own line as a single subfeature. These would be
linked together by a shared attribute (`transcript_id` in the case of `GTF2`_;
`parent` in the case of `GFF3`_) to reconstruct the parental transcript.

A parser cannot know whether it has collected all of the sub-features
needed to assemble a discontinuous feature until it receives information
determining this is so. This information could be:

  - In a `GFF3`_ file, the special line::
    
        # this line is a comment, ignored by GFF3 parsers.
        ###
        # the line above is not a comment, but a GFF3 instruction!
        # this line and the line above it are comments. 
        
    which indicates all features in memory may be assembled.
  - In a sorted `GTF2`_ or `GFF3`_ file, a change in chromosomes, indicating
    all features on the previous chromosome may be assembled.
  - The end of the annotation file 

In all cases, the parser has to hold all collected features in memory until
it it receives some signal that it is ok to assemble the features. This costs
memory, time, and disk space, but allows subfeatures to have their own
annotation data, and arbitrary keywords.

However, if all of your features are continuous, they can all represented one
a single line, and don't need to be assembled. In this case, a `GTF2`_ or `GFF3`_
poses little additional cost.


In summary
..........
The table below summarizes the discussion above: 

==========   =====================================    ==========================    ======================   ==============
**Format**   **Features that are not transcripts**    **Multiple feature types**    **Feature attributes**   **Memory use**
             **or parts of transcripts**    
----------   -------------------------------------    --------------------------    ----------------------   --------------
`BED`_       Yes                                      No                            No                       Low

`BED+X`_     Yes                                      If specified in extra         1 per extra column       Low
                                                      column
                                                      
`BigBed`_    Yes                                      If specified in extra         1 per extra column       Low
                                                      column
                                                      
`GTF2`_      No                                       Yes                           Unlimited                High for discontinuous features

`GFF3`_      Yes                                      Yes                           Unlimited                High for discontinuous features
==========   =====================================    ==========================    ======================   ==============


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
  - The `UCSC file format FAQ`_ for details on file formats and further discussion
    of their capabilities, advantages, and disadvantages
  - The `GFF3 specification <GFF3>`_ for details on GFF3 files
  - :doc:`/concepts/coordinates` for information on genomic coordinates

