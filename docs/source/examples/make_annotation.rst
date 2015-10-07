Creating custom `BED`_, `GTF2`_, and `GFF3`_ files
==================================================

In this tutorial, we describe how to make custom :term:`genome annotations <annotation>`
in `BED`_ and `GTF2`_ formats. These can then be used as any other annotation file:

  - to view custom regions of interest in a :term:`genome browser` like `IGV`_ or `UCSC`_
  - to perform analyses -- using pipelines in :data:`plastid` or other programs -- on
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
    
 #. Use the :meth:`~plastid.genomics.roitools.SegmentChain.as_bed`,
    :meth:`~plastid.genomics.roitools.SegmentChain.as_gtf`, and 
    :meth:`~plastid.genomics.roitools.SegmentChain.as_gff3` methods to export
    the |SegmentChains| or |Transcripts| in the desired format. To make
    `BigBed`_ files, see :ref:`make-annotation-bigbed`.

:data:`plastid` automatically handles all coordinate conversions, splicing,
text formatting, and attribute conversion.


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

     #. Combine the sequences of the human genome and the reporter construct
        (lenti223.fa) into a single `FASTA`_ file, and build an alignment
        index for `bowtie`_. From the terminal:

         .. code-block:: shell

            $ cat hg38Patch3.fa lenti223.fa >combined_sequences.fa
            $ bowtie-build combined_sequences.fa my_combined_index

     #. Create a custom annotation describing the coordinates of the reporter gene
        with respect to the vector sequence. Coordinates should be :term:`0-indexed`
        and :term:`half-open` (i.e. typical Python idioms). This can be done 
        interactively in Python by creating |Transcripts| describing the reporter::

            >>> from plastid.genomics.roitools import GenomicSegment, SegmentChain, Transcript

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
DNA-binding protein dCAS9 (see :cite:`Gilbert2014`).
This requires us to:

  #. Define windows upstream of the transcription start sites (TSS) for the genes we
     wish to knock down. This step we'll perform here.

  #. Search within those windows for genomic sequences that we can target with guide
     RNAs for dCAS9. For details on that procedure, see :cite:`Gilbert2014`.

We'll use the same transcript annotation as in the example above, but first we'll download the
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
    >>> from plastid.readers.gff import GTF2_TranscriptAssembler
    >>>
    >>> transcripts = GTF2_TranscriptAssembler(open("gencode.v23.annotation.gtf"),sorted=True)
    >>>
    >>> for tx in transcripts:
    >>>     chrom, tx_start, strand =  tx.get_genomic_coordinate(0)
    >>>     # for plus-strand transcripts, TSS is 5' of transcript on chromosome
    >>>     if strand == "+":
    >>>         tss_window = GenomicSegment(chrom,tx_start-500,tx_start,strand)
    >>>     # for minus-strand transcripts, TSS is 3' of transcript on chromosome
    >>>     elif strand == "-":
    >>>         tss_window = GenomicSegment(chrom,tx_start,tx_start+500,strand)
    >>>
    >>>     # write coordinates of TSS to a BED file
    >>>     tss_window = SegmentChain(tss_window,ID="%s_tss_window" % tx.get_name())
    >>>     bed_out.write(tss_window.as_bed())
    >>>
    >>>     # write genomic sequence of TSS to a FASTA file
    >>>     tss_window_sequence = tss_window.get_fasta(genome)
    >>>     seq_out.write(tss_window_sequence)
    
    >>> # close files
    >>> bed_out.close()
    >>> seq_out.close()

The `fasta`_ file of sequences can then be processed with any pipeline, and the
TSS windows viewed in a :term:`genome browser`, like `IGV`_ or the `UCSC genome browser`_.



 .. _make-annotation-bed-xplusy:

Making custom :term:`BED X+Y` files
-----------------------------------

To export attributes of a |SegmentChain| or |Transcript| as extra columns
in a :term:`BED 12+Y` format, pass the `extra_columns` keyword to the
:meth:`plastid.genomics.roitools.SegmentChain.as_bed` method::

    >>> attr = { "ID" : "some feature ID",
                 "extra_field_1" : 542,
                 "extra_field_2" : "some extra field",
               }

    >>> my_chain = Transcript(GenomicSegment("chrA",100,150,"+"),
                              GenomicSegment("chrA",500,550,"+"),
                              **attr)
    >>> my_chain.as_bed(extra_columns=["extra_field_1","extra_field_2"])
    chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,

    >>> my_chain.as_bed(extra_columns=["extra_field_1","extra_field_2"])
    chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,    542    some extra field

If an attribute is not defined, the column will be left empty ""::

    >>> my_chain.as_bed(extra_columns=["extra_field_1","nonexistent_field","extra_field_2"])
    chrA    100    550    some feature ID    0    +    100    100    0,0,0    2    50,50,    0,400,    542        some extra field





 .. _make-annotation-bigbed:

Making `BigBed`_ files
----------------------
`BigBed`_ files are easily made from `BED`_ files using `Jim Kent's utilities`_.
To make a `BigBed`_ file:

 #. Create a custom `BED`_ or :term:`BED X+Y`, file, following the examples above.
    For :term:`BED X+Y` files, consider making an optional `autoSql`_ description
    of the names & data types of the extra columns. This will allow parsers to 
    convert these to native types when reading the `BigBed`_ file.
 
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


-------------------------------------------------------------------------------


See also
--------
  - :ref:`data-annotation-format` for a brief overview of the costs & benefits
    of `BED`_, `BigBed`_, `GTF2`_ and `GFF3`_ files.
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
  - `Jim Kent's utilities`_ for more info on making `BigBed`_ files.


