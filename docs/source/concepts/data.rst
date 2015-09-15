Categories and formats of genomics data
=======================================
This document contains very, very brief overviews of the types of data
encountered in genomics, and of some common file formats.

.. _quickstart-data:

Types of data
-------------

Generally speaking, genomics data comes in a few flavors:

    Sequence
        The nucleotide sequence of a chromosome, contig, transcript,
        or a set of these. These are typically maintained by public databases,
        such as `UCSC <UCSC genome browser>`_, `Ensembl`_, and `RefSeq`_. The
        genome sequence for a given organism frequently is available in several
        editions, called :term:`builds <genome build>` or :term:`assemblies <genome build>`.
    
    :term:`Annotations <annotation>`
        Descriptions of features -- e.g. genes, transcripts, SNPs, start codons
        -- that appear in genomes or transcripts. Annotations typically include
        coordinates (chromosome name, chromosome positions, and a chromosome
        strand), as well as properties (gene name, function, GO terms, et c) of
        a given feature.
        
        :term:`Annotations <annotation>` are maintained by the same public
        databases that maintain sequence information, because the coordinates
        in each annotation are specific to the :term:`genome build` upon which
        it is based. In other words, annotations and sequences must be matched!
        Pay particular attention to this when downloading annotations as 
        supplemental data from a journal article.
        
    Quantitative data
        Any kind of numerical value associated with a chromosomal
        position. For example, the degree of phylogenetic conservation between a 
        set of organisms, at each position in the genome. Or, the strength of 
        transcription factor binding to a chromosomal position in a ChIP-seq dataset.
        
        Because quantitative data associates values with chromosomal coordinates,
        it can be considered an :term:`annotation` of sorts. It is therefore
        important again to make sure that the coordinates in your data file
        match the :term:`genome build` used by your feature :term:`annotation`
        and/or :term:`read alignments`.
        
    :term:`Read alignments <read alignments>`
        A record matching a short sequence of DNA to a region of identical or similar
        sequence in a genome. In a :term:`high-throughput sequencing` experiment,
        alignment of short reads identifies the genomic coordinates from which
        each read presumably derived.
        
        :term:`Read alignments <read alignments>` can be produced by running
        sequencing data through alignment programs,
        such as `Bowtie`_, `Tophat`_, or `BWA`_. 
        
        :term:`Read alignments <read alignments>`
        can be converted to quantitative data by applying a :term:`mapping rule`,
        that converts the reads to :term:`counts`. For example, one could count
        the number of 5' ends of reads that align to each position in a genome. For
        an in-depth discussion of :term:`mapping rules <mapping rule>`, see
        :doc:`concepts/mapping_rules`.


Formats of data
---------------
One of the design goals of :data:`plastid` is to insulate users from the esoterica
of the various file formats used in genomics. But, two points are relevant:

  #. It is important for users to recognize the file types names in order to 
     identify the files they have or need to download.
     
  #. Some file formats are *indexed* and others are not. Indexed files are
     memory-efficient, because computer programs don't need to read the entire
     file to find the data of interest; instead, they can read the index and
     just fetch the desired portion of the data.
     
     However, indexed files are frequently compressed, which can make reading them 
     slower to parse. For small genomes that don't use much memory in the first
     place (e.g. yeast, *E. coli*), the meager memory savings aren't worth this
     speed cost. The exception is for short :term:`read alignments`, where indexed
     `BAM`_ files are universally recommended. 

.. TODO later: update when format support changes

Below is a table of commonly used file formats. At present, :data:`plastid` handles
all of these except `BigWig`_, either natively or via `Pysam`_ (`BAM`_ files),
`Biopython`_ (`FASTA`_), or `2bitreader`_ (`2bit <twobit>`_).

    =====================   ===================================   ===================
    **Data type**           **Unindexed formats**                 **Indexed formats**
    ---------------------   -----------------------------------   -------------------
    Sequence                `FASTA`_                              `2bit <twobit>`_
    
    Annotations             `BED`_, `GTF2`_, `GFF3`_, `PSL`_      `BigBed`_ 
    
    Quantitative data       `bedGraph`_, `wiggle`_                `BigWig`_
    
    Read alignments         `bowtie`_, `SAM`_, `PSL`_             `BAM`_ 
    =====================   ===================================   ===================
 
 
`BED`_, `GTF2`_, `GFF3`_, and `PSL`_ files can be indexed via `tabix`_.
:data:`plastid` supports (via `pysam`_) reading of `tabix`_-compressed files too.


Why are there so many formats?
------------------------------

There are a number of answers to this:

 #. Genomics is a young science, and for a long time there was no consensus
    on how best to store data. This dialogue is, in fact, still ongoing.
     
 #. It became apparent that file formats that work well with small genomes
    become very onerous for mammalian-sized genomes. This is why, for example,
    the `2bit <twobit>`_, `BigBed`_, and `BigWig`_ formats were created. 

 #. The various file formats have their own strengths and weaknesses. These
    are detailed in :ref:`data-annotation-format`
    

 .. _data-annotation-format:

Which annotation format should I use?
-------------------------------------
When choosing a feature annotation format, consider the following questions:

  - Will the annotation contain features that are not transcripts?
  - Will multiple types of features be stored in the same file?
  - Does rich attribute information need to be saved in the file?
  - Are features discontinuous?
  - Is the computing environment limited for processing power or memory and/or
    is the feature annotation very large?


`BED`_, :term:`BED X+Y`, & `BigBed`_
....................................
`BED`_-family files contain a single record per line. And, in contrast
to `GTF2`_ or `GFF3`_ files, single records -- like transcripts -- can
be discontinuous. This makes `BED`_ files computationally
cheap to parse, because each line is a complete record. In contrast, 
in `GTF2`_ and `GFF3`_ files, discontinuous features like transcripts need
to be assembled from multiple continuous records (e.g. records describing
individual exons).

`BED`_ files contain columns that describe only the following attributes:

  - feature name
  - feature coordinates (feature can be discontinuous, like a multi-exon transcript)
  - feature coding region start & stop  
  - a score for the feature
  - a color for rendering the feature in a genome browser

Note that *there is no attribute for feature type:* typically all records
in a `BED`_ file are of the same type (e.g. every record is a transcript
or an alignment or a ChIP binding site, et c).

`BigBed`_ and :term:`BED X+Y` formats can include additional attributes in additional
columns, but every entry in each column must be the same type of attribute 
(e.g. a "gene id" column can only contain gene IDs).



`GTF2`_ & `GFF3`_
..................
Unlike `BED`_-family files, `GTF2`_ and `GFF3`_ files are hierarchical:
features have parents and children, which themselves are other features. Continuous
features are represented on a single line. Discontinuous features -- like
transcripts -- are represented on multiple lines -- for example, one
line per exon, one line per intron, and one line per continous portion
of a coding region. These sub-features are linked together via parent-child
attributes (`'Parent'` in for `GFF3`_; `'gene_id'` and `'transcript_id'` in
`GTF2`_), which associate them with the discontinuous feature they represent.

This has several important implications:

 #. Sub-features in `GTF2`_ & `GFF3`_ can have their own attributes,
    which differ from the attributes of their parent features.
 
 #. In order to reconstruct a discontinuous feature like a transcript,
    `GTF2`_ & `GFF3`_ parsers need to collect all of the required subfeatures.
    However, parsers only know when they have collected all of the required features
    if they receive information indicating this is so. This information could be:

      - In a `GFF3`_ file, the special line `'###'`::
        
            ###
            # the line above is not a comment, but a GFF3 instruction!
            # this line and the line above it are comments.  
            
        which indicates all features in memory may be assembled.
      - In a sorted `GTF2`_ or `GFF3`_ file, a change in chromosomes, indicating
        all features on the previous chromosome may safely be assembled.
      - The end of the annotation file 

    In all cases, a `GTF2`_ or `GFF3`_ parser has to hold a potentially large
    set of subfeatures in memory until it it receives some signal that all related
    subfeatures have been collected. This costs memory, time, and disk space, and
    can become unwieldy for large genomes.

However,
a major advantage of `GTF2`_ and `GFF3`_ files is that they contain a column (column 9)
for arbitrary key-value pairs of attributes (such as GO terms, descriptive paragraphs,
IDs that cross-reference different databases). This allows different features to have
different types of attributes.

The primary difference between `GTF2`_ and `GFF3`_ formats is that, formally, 
`GTF2`_ files only describe transcripts and their parts, according to a defined
schema. The complete list of valid record types in `GTF2`_ is:

  - CDS
  - start_codon
  - stop_codon
  - 5UTR
  - 3UTR
  - inter
  - inter_CNS
  - intron_CNS
  - exon

`GFF3`_ files can describe any type of feature with any schema of parent-child
hierarchy. This makes `GFF3`_ the most flexible format. The cost of this
flexibility is that, without knowing the parent-child schema, `GFF3`_ parsers 
don't know which chidl subfeatures to assemble into complex parent features.

Similarly, given an assembled feature in Python (represented as a
|SegmentChain|), in the absence of a schema there is ambiguity surrounding
what types the parent |SegmentChain| and each of its children (|GenomicSegments|)
should be rendered as in `GFF3`_ output. Due to this ambiguity, attempts to call the
:meth:`~plastid.genomics.roitools.SegmentChain.as_gff3` method on a multi-segment
|SegmentChain| will raise an :py:obj:`AttributeError`.

 .. _data-export-gff3:

Instead, users may export the individual features from which the
multi-segment |SegmentChain| was constructed, setting `'ID'`, `'Parent'`,
and `'type'` attributes in each child feature's `attr` dict::

    >>> # a multi-segment chain
    >>> my_alignment
    <SegmentChain segments=2 bounds=chrI:212353-214802(+) name=some_alignment>
    >>> my_alignment.attr
    {'ID': 'some_alignment', 'type': 'alignment'}
    >>> list(my_alignment)
    [<GenomicSegment chrI:212353-212900 strand='+'>,
     <GenomicSegment chrI:214313-214802 strand='+'>]
    >>> my_alignment.as_gff3()
    # AttributeError!

    >>> # make a single, continuous feature with the endpoints of `my_alignment`
    >>> # 'ID' attribute should match 'ID' of my_alignment
    >>> alignment_span = SegmentChain(my_alignment.spanning_segment,ID="some_alignment",type="alignment")

    >>> # then make a subfeature for each segment `my_alignment`,
    >>> # 'Parent' attribute should match the 'ID' attribute of `alignment_span`
    >>> block1 = SegmentChain(my_alignment[0],Parent=my_alignment.get_name(),type="aligned_block")
    >>> block2 = SegmentChain(my_alignment[1],Parent=my_alignment.get_name(),type="aligned_block")

    >>> # write to file
    >>> features = [alignment_span,block1,block2]
    >>> with open("some_file.gff","w") as gff_out:
    >>>     for feature in features:
    >>>         gff_out.write(feature.as_gff3())

In contrast, multi-segment |Transcripts| *can* be unambiguously exported to `GFF3`_;
they are rendered using the ontology from 
`Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_.


In summary
..........
The table below summarizes the discussion above: 

===============   =====================================    ==========================    ======================   ==============
**Format**        **Features that are not transcripts**    **Multiple feature types**    **Feature attributes**   **Memory use**
                  **or parts of transcripts**    
---------------   -------------------------------------    --------------------------    ----------------------   --------------
`BED`_            Yes                                      No                            No                       Low

:term:`BED X+Y`   Yes                                      If specified in extra         1 per extra column       Low
                                                           column
                                                      
`BigBed`_         Yes                                      If specified in extra         1 per extra column       Low (and indexed)
                                                           column
                                                      
`GTF2`_           No                                       Yes                           Unlimited                High for discontinuous features

`GFF3`_           Yes                                      Yes                           Unlimited                High for discontinuous features
===============   =====================================    ==========================    ======================   ==============


Getting the most out of your time & data
----------------------------------------

Starting a new type of analysis is rarely straightfoward. But, it is possible 
to save some time by following several practices:

 #. Make sure your :term:`annotation` matches your :term:`genome build`. e.g.
    do not use the *mm9* mouse genome annotation with the *mm10* sequence
    assembly. Do not mix `Ensembl`_'s human genome build *GRCh38* and
    `UCSC <UCSC genome browser>`_'s similar-but-still-different *hg38*.

 #. If using a large genome (e.g. *Drosophila* or larger), consider using
    non-hierarchical (e.g. `BED`_) and possibly indexed (e.g. `BigBed`_,
    `BigWig`_ ) formats instead of non-indexed formats.

 #. Work from alignments in `BAM`_, rather than `bowtie`_, format.

-------------------------------------------------------------------------------

See also
--------
  - :class:`~plastid.genomics.roitools.SegmentChain` and
    :class:`~plastid.genomics.roitools.Transcript` for details on these classes
  - The `UCSC file format FAQ`_ for details on file formats and further discussion
    of their capabilities, advantages, and disadvantages
  - The `GFF3 specification <GFF3>`_ for details on GFF3 files
  - :doc:`/concepts/coordinates` for information on genomic coordinates
  - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_,
    for a description of a common `GFF3`_ feature ontology
  - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_,
    for the current SO consortium release.


