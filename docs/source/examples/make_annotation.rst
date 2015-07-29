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
made manually in a text editor or spreadsheet. However, an easier and
potentially less error-prone process is to:

 #. Enter an interactive Python session
 
 #. Manually construct |SegmentChains| and/or |Transcripts| that represent
    the :term:`features <feature>` of interest, using :term:`0-indexed`,
    :term:`half-open` coordinates
    
 #. Use the :meth:`~yeti.genomics.roitools.SegmentChain.as_bed`,
    :meth:`~yeti.genomics.roitools.SegmentChain.as_gtf`, and 
    :meth:`~yeti.genomics.roitools.SegmentChain.as_gff3` methods to export
    the |SegmentChains| or |Transcripts| in the desired format. To make
    `BigBed`_ files, see :ref:`make-annotation-bigbed`.

All coordinate conversions, splicing, text formatting, and attribute
conversion is handled for you.


Examples
--------

A Reporter construct
....................
Suppose we made a cell line containing a reporter construct, and that we 
want to measure the expression level of this reporter along with the expression
levels of all other genes in an :term:`RNA-seq` experiment. In this example:

  - The host cell line is K562
  - Reporters are introduced by a lentiviral vector called 'lenti223'
  - The vector contains two reporters, both on the plus-strand of the virus:
    
      - Superfolder GFP driven by some promoter
      - mCherry driven by some other promoter


An analysis pipeline for this experiment might look the following:

     #. Obtain a genome :term:`annotation` for the host genome. K562 cells
        are human, so we'll obtain the human *hg38p3* reference sequence 
        from `UCSC`_, and the matching annotation from `GENCODE`_.
        From the terminal:

         .. code-block:: shell

            # get chromosomal sequence in FASTA format for making bowtie index
            $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/hg38Patch3/hg38Patch3.fa.gz

            # download annotation
            $ wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz

            # unzip them
            $ gunzip *gz

     #. Combine the sequences into a single file. From the terminal:

         .. code-block:: shell

            $ cat hg38Patch3.fa lenti223.fa >combined_sequences.fa
            $ bowtie-build combined_sequences.fa my_combined_index

     #. Create a custom annotation describing the coordinates of the reporter gene
        with respect to the vector sequence. Coordinates should be :term:`0-indexed`
        and :term:`half-open` (e.g. typical Python idioms). We do this in an
        interactive Python session, by creating |Transcripts| describing the reporter::

            >>> from yeti.genomics.roitools import GenomicSegment, SegmentChain, Transcript

            # GFP transcript, containing 100 bp of 5' UTR and 150 bp of 3' UTR
            # 714bp coding region from bases 945-1659
            >>> gfp = Transcript(GenomicSegment("lenti223",845,1809,"+"),ID="sfGFP",cds_genome_start=945,cds_genome_end=1659)

            # mCherry transcript, similarly constructed
            >>> rfp = Transcript(GenomicSegment("lenti223",2100,3061,"+"),ID="mCherry",cds_genome_start=2200,cds_genome_end=2911)

            # now, write out features 

            >>> with open("custom.gtf","w") as fout:
            >>>     fout.write(gfp.as_gtf())
            >>>     fout.write(rfp.as_gtf())
            >>>     fout.close()

        The file ``custom.gtf`` should look something like this:

         .. code-block:: shell

            lenti223    .    exon           846     1809    .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; ID "sfGFP";
            lenti223    .    CDS            946     1656    .    +    0    gene_id "gene_sfGFP"; transcript_id "sfGFP"; ID "sfGFP";
            lenti223    .    start_codon    946     948     .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; cds_start "100"; cds_end "814"; ID "sfGFP";
            lenti223    .    stop_codon     1657    1659    .    +    .    gene_id "gene_sfGFP"; transcript_id "sfGFP"; cds_start "100"; cds_end "814"; ID "sfGFP";
            lenti223    .    exon           2101    3061    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; ID "mCherry";
            lenti223    .    CDS            2201    2908    .    +    0    gene_id "gene_mCherry"; transcript_id "mCherry"; ID "mCherry";
            lenti223    .    start_codon    2201    2203    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; cds_start "100"; cds_end "811"; ID "mCherry";
            lenti223    .    stop_codon     2909    2911    .    +    .    gene_id "gene_mCherry"; transcript_id "mCherry"; cds_start "100"; cds_end "811"; ID "mCherry";


        Then, merge the annotations, from the terminal:

         .. code-block:: shell

            $ cat gencode.v23.annotation.gtf custom.gtf >my_cell_line_combined.gtf

     #. Align data from ``some_file.fq`` in `tophat`_. From the terminal:

         .. code-block:: shell

            $ tophat -G my_cell_line_combined.gtf -o my_alignments my_combined_index some_file.fq


     #. Perform quantitation using pipeline of choice (e.g. `cufflinks`_, |cs|, |counts_in_region|, or something else)
        


Identifying target sites for gene knockdown via dCAS9
.....................................................
Suppose we wish to knock down target genes using the programmable
DNA-binding protein dCAS9 protein (dCAS9; see :cite:`Gilbert2014`).
This requires us to:

  #. Define windows upstream of the transcription start sites (TSS) for the genes we
     wish to knock down. This step we'll perform here.

  #. Search within those windows for genomic sequences that we can target with guide
     RNAs for dCAS9. For details on how to do that, see :cite:`Gilbert2014`. Here,
     we'll just fetch the nucleotide sequence.

It is possible to download known annotations of transcription start sites. For the sake of this
example, let's suppose those files weren't available, and that we'd need to define these on our own.

We'll use the same annotation as in the example above, but first we'll download the
`2bit <twobit>`_-formatted version of the genome sequence, which requires less memory to read
(*n.b.* if you haven't already, you need to install the
`twobitreader <https://pypi.python.org/pypi/twobitreader>`_ package from `PyPI`_).
From the terminal:

 .. code-block:: shell

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/hg38Patch3/hg38Patch3.2bit


Then, within a Python session, read each transcript and create TSS windows::

    >>> # open genome sequence
    >>> from twobitreader import TwoBitFile
    >>> genome = TwoBitFile("hg38Patch3.2bit")

    >>> # open files where we'll save our data
    >>> bed_out = open("tss_windows.bed","w")
    >>> seq_out = open("tss_window_sequences.fa","w")

    >>>  # open transcripts and process one-by-one
    >>> from yeti.readers.gff import GTF2_TranscriptAssembler
    >>> transcripts = GTF2_TranscriptAssembler(open("gencode.v23.annotation.gtf"),sorted=True)
    >>> for tx in transcripts:
    >>>     chrom, tx_start, strand =  tx.get_genomic_coordinate(0)
    >>>     if strand == "+":
    >>>         tss_window = GenomicSegment(chrom,tx_start-500,tx_start,strand)
    >>>     elif strand == "-":
    >>>         tss_window = GenomicSegment(chrom,tx_start,tx_start+500,strand)
    >>>
    >>>     tss_window = SegmentChain(tss_window,ID="%s_tss_window" % tx.get_name())
    >>>     bed_out.write(tss_window.as_bed())
    >>>
    >>>     tss_window_sequence = tss_window.get_fasta(genome)
    >>>     seq_out.write(tss_window_sequence)
    
    >>> bed_out.close()
    >>> seq_out.close()

The `fasta`_ file of sequences can then be processed with any pipeline, and the
TSS windows viewed in a :term:`genome browser`, like `IGV`_ or the `UCSC genome browser`_.


 .. _make-annotation-choose-format:

Choosing a format: `BED`_, `BigBed`_, `GTF2`_, and `GFF3`_
----------------------------------------------------------

Which format to choose depends on your purposes. The following questions
highlight the differences between formats:

Are your features something other than exons, coding regions, UTRs or transcripts?
..................................................................................
If so, you need to use something other than `GTF2`_, which, according to its
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
If so, use something other than `BED`_ format. In classic `BED`_ formats,
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


 .. Note::

    The `GFF3`_ specification allows any schema of parent-child hierarchy,
    making `GFF3`_ files incredibly flexible.
    
    However, |SegmentChains| are unaware of which schema is in use at any given moment,
    and therefore do not know what types the parent |SegmentChain| and each
    of its children (|GenomicSegments|) should be rendered as in `GFF3`_ output.
    Due to this ambiguity, attempts to call the :meth:`~yeti.genomics.roitools.SegmentChain.as_gff3`
    method on a |SegmentChain| that requires parent-child relationships for export
    -- i.e. all multi-segment chains -- will raise an :py:obj:`AttributeError`.
    
    Instead, users may export the individual features from which the
    multi-segment |SegmentChain| was constructed, setting *ID*, *Parent*,
    and *type* attributes in each child feature's `attr` dict::

        >>> # a multi-segment chain
        >>> my_alignment
        <SegmentChain segments=2 bounds=chrI:212353-214802(+) name=some_alignment>
        >>> my_alignment.attr
        {'ID': 'some_alignment', 'type': 'alignment'}
        >>> list(my_alignment)
        [<GenomicSegment chrI:212353-212900 strand='+'>,
         <GenomicSegment chrI:214313-214802 strand='+'>]

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


Does your dataset include multiple types of features?
.....................................................
If so, use `GTF2`_ or `GFF3`_. Because `BED`_ files contain no column to
describe feature type, it is simplest to make sure all features in the `BED`_
file are of a single type.


Are the features you care about discontinuous? And is your computer limited for memory?
.......................................................................................
If so, use one of the `BED`_-family formats. In `BED`_ files, each feature
-- even discontinuous
features like multi-exon transcripts -- are represented as single lines. This
means that programs don't need to search through a file to find all of the
pieces (e.g. exons, piecse of coding regions, et c) that make up a feature.

In contrast, in `GTF2`_ and `GFF3`_ files, each line can only contain a continuous
feature or sub-feature. So, to represent a multi-exon transcript, each exon
would be represented on its own line as a single subfeature. These would be
linked together by a shared attribute (`'transcript_id'` in the case of `GTF2`_;
`'parent'` in the case of `GFF3`_) to reconstruct the parental transcript.

A `GTF2`_ or `GFF3`_ parser cannot know whether it has collected all of the sub-features
needed to assemble a discontinuous feature until it receives information
indicating this is so. This information could be:

  - In a `GFF3`_ file, the special line::
    
        # this line is a comment, ignored by GFF3 parsers.
        ###
        # the line above is not a comment, but a GFF3 instruction!
        # this line and the line above it are comments. 
        
    which indicates all features in memory may be assembled.
  - In a sorted `GTF2`_ or `GFF3`_ file, a change in chromosomes, indicating
    all features on the previous chromosome may be assembled.
  - The end of the annotation file 

In all cases, a `GTF2`_ or `GFF3`_ parser has to hold all collected features in memory until
it it receives some signal that all related features have been collected. This costs
memory, time, and disk space, but allows subfeatures to have their own
annotation data, and arbitrary keywords.

**However**, if all of your features are continuous, they can all represented one
a single line in `GTF2`_ and `GFF3`_, and don't need to be assembled. In this case,
`GTF2`_ or `GFF3`_ formats pose no additional cost compared to `BED`_.


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

 .. _make-annotation-bigbed:

Making `BigBed`_ files
----------------------
`BigBed` files are easily made from `BED`_ files using `Jim Kent's utilities`_.
To make a `BigBed`_ file:

 #. Create a custom `BED`_ or `BED+X`_, file, following the examples above

 #. Sort the `BED`_ file by chromosome and start position. This is easily 
    done in a terminal session:
    
     .. code-block:: shell

        $ sort -k1,1n -k2,2n my_annotation.bed >my_annotation_sorted.bed

 #. Download and install `Jim Kent's utilities`_, which include the
    ``bedToBigBed`` program.

 #. Obtain a chromosome/contig ``.sizes`` file. If using genome builds from
    `UCSC`_, these can be downloaded using the ``fetchChromSizes`` program
    included with `Jim Kent's utilities`_. For example:

     .. code-block:: shell

        $ fetchChromSizes hg38 >>hg38.sizes 

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
  - `Sequence Ontology (SO) v2.53 <http://www.sequenceontology.org/browser/>`_,
    for a description of a common `GFF3`_ feature ontology
  - `SO releases <http://sourceforge.net/projects/song/files/SO_Feature_Annotation/>`_,
    for the current SO consortium release.


