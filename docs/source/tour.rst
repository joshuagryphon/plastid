Getting started with ``yeti``
====================================

What you need
-------------
 	* Alignments
 	* Annotations


Command-line scripts
--------------------
``yeti`` implements a number of common sequencing workflows end-to-end,
including several workflows specific to analysis of ribosome profiling.
These include, among others:

    * :doc:`Measuring read density </generated/yeti.bin.cs3>` in regions
      of interest, for example, to obtain gene expression values
    
    * :doc:`Determining P-site offsets </generated/yeti.bin.psite>` for
      ribosome read mapping
    
    * :doc:`Performing metagene analyses </generated/yeti.bin.metagene>` for
      any sort of sequencing data
    
    * :doc:`Creating browser tracks </generated/yeti.bin.make_wiggle>` 
      from BAM or bowtie files, after applying optional read mapping functions
      (e.g. for P-site assignment, in the case of ribosome profiling data) 
    
    * :doc:`Measuring sub-codon phasing </generated/yeti.bin.phase_by_size>`
      for ribosome profiling data


For a complete list, as well as examples, see the section on :doc:`Command-line scripts </cli_howto>`


Walk-throughs
-------------
``yeti`` is useful for both command-line and interactive analyses of
sequencing data. We have included a number of :doc:`How-tos </how_to_list>`
to help you get started, including:

	* how-to number 1
	
	* how-to number 2

	* how-to number 3



.. _overview-of-data-structures:

Overview of data structures
---------------------------
``yeti`` includes a number of data structures that facilitate ad-hoc
analyses. 

.. _ivcollection_description :

|GenomicSegment|
^^^^^^^^^^^^^^^^^
	A |GenomicSegment| models a single, continuous region of interest of the genome,
	and is the fundamental building block of more complex features. 
	A |GenomicSegment| can be as large as a chromosome, or as small as a nucleotide. 
	|GenomicSegment| s provide convenience methods for evaluating overlap
	or containment with other features.

|SegmentChain| & |Transcript|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	|SegmentChain| & its subclass |Transcript| model discontinuous features of
	the genome, such as gapped alignments, chromosomal breaks, splice junctions,
	or transcripts. |SegmentChain| and its subclasses provide methods for:
	
	    - converting between |SegmentChain|-centric and genomic
	      coordinate spaces
	    - fetching vectors of values (e.g. high-throughput sequencing counts or 
	      conservation data) at each position in the collection, moving from
	      the fiveprime toward the threeprime end and accounting
	      for splicing if the SegmentChain is discontinuous
	    - testing equality to, overlap with, or containment of other |SegmentChain| s
	      or |GenomicSegment| s et c
	    - importing/exporting to/from various annotation formats 
	      (e.g. BED, BigBed, GTF, GFF, or PSL), using appropriate
	      :doc:`readers </generated/yeti.readers>`
	     
	     
	For example, to load transcripts from a BED file::
	 
	    >>> from yeti.readers.bed import BED_to_Transcripts
	    >>> my_transcripts = BED_to_Transcripts(open("some_file.bed"))
	    >>> for transcript in my_transripts:
	            # do_something
	
	
	Fetch the fiveprime 200 nucleotides of a transcript, regardless of whether
	it is on the plus or minus-strand::
	
	    >>> new_ivc = my_transcript.get_subchain(0,200)
	
	
	Similarly, we can ask what a genomic coordinate is, relative to an
	|SegmentChain|::
	
	    >>> transcript_coordinate = my_transcript.get_segmentchain_coordinate(chrom,genomic_coordinate,strand)
	
	   
	or to fetch the coding region of a transcript (if it is coding)::
	
	    >>> my_transcript.get_cds()
	        # OUTPUT HERE
	
	
	|SegmentChain| and its subclasses can also fetch their own sequences from dictionaries
	(or dictionary-like objects). These sequences will automatically be spliced if
	the |SegmentChain| has several exons or sub-regions, and reverse-complemented
	as necessary::
	
	    >>> my_transcript.get_sequence(dict_of_chrom_sequences)
	        "tcgataccatacgtgcactgaagarta"
	
	
	They can also fetch vectors of sequencing counts from objects called |GenomeArray| s,
	again accounting for splicing, so that each position in the returned vector corresponds to a position
	in the |SegmentChain|, from the fiveprime to the threeprime
	end::
	
	    >>> my_transcript.get_counts(genome_array)
	        [3,5,1,4,6, ... ]
	
	
	Fore more information, see the documentation for |SegmentChain|,
	|Transcript|, and the :py:mod:`~yeti.genomics.roitools` module.
        
 
|GenomeArray|, |BAMGenomeArray|, & related subclasses
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    |GenomeArray| s store count data at each position in the genome. Data can be
    imported from count files (e.g. `Wiggle`_, `bedGraph`_) as well as alignment files
    (in `bowtie`_ or `BAM`_ format). For very large genomes a sparse implementation
    is provided by |SparseGenomeArray|. 
    
    When importing alignment files, users can specify arbitrary :doc:`mapping functions </mapping_functions>`
    that determine how reads should be converted into counts (e.g., to their
    fiveprime ends, threeprime ends, or, somewhere in between, optionally
    as a function of read length). ``yeti`` already includes mapping
    functions to map read alignments:
    
        * to their fiveprime or threeprime ends, with or without offsets from
          the end of the read. These offests can be constant, or a function of 
          read length (as is typical for ribosome profiling data). 
         
        * fractionally over their entire lengths (e.g. for RNA-seq)
        
        * fractionally to all positions covered by a central portion of the read
          alignment, after excluding a user-defined number of positions on each
          send of the read (as in ribosome profiling data from *E. coli* or *D. melanogaster)*.


    For further information, see:
    
		* The discussion of :doc:`mapping functions </mapping_functions>`
    	* The module documentation for :py:mod:`~yeti.genomics.genome_array`
    	* Class documentation for |GenomeArray|, |BAMGenomeArray|, and |SparseGenomeArray|


|GenomeHash| and |BigBedGenomeHash|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    It is frequently useful to retrieve features that overlap specific regions 
    of interest in the genome, for example, to find transcripts that overlap one
    another. However, it would be inefficient to compare features that are
    too far apart in the genome to overlap in the first place. 
    
    |GenomeHash| and |BigBedGenomeHash| index genomic annotations by location
    to avoid making unnecessary comparisons. A |GenomeHash| may be created  
    from a list or dictionary of features (e.g. |SegmentChain| s or |Transcript| s),
    or directly loaded from a genome annotation (in BED, GTF2, GFF3, or PSL format).
    
    A |BigBedGenomeHash| may be created from a BigBed file, and takes advantage
    of the indices already present in the BigBed file to avoid loaded annotations
    into memory before they are used (if they even are at all).
     
     
    For example, to find all features between bases 10000-20000 on the plus strand of
    chromosome *chrI*::
    
        >>> my_hash = GenomeHash(list_or_dict_of_transcripts)
        >>> roi = GenomicSegment("chrI",10000,20000,"+")
        >>> features_overlapping_my_roi = my_hash[roi]
    
    
    Or, to find features that overlap a |Transcript| or |SegmentChain|::
    
        >>> overlapping_transcripts = my_hash[my_transcript]
    
    
    For more information, see the module documentation for 
    :py:mod:`~yeti.genomics.genome_hash`.



See also
--------
For more documentation, please see:

	* :ref:`API documentation <modindex>`
	
	* :doc:`Complete list of command-line scripts </cli_howto>`
	
	* :doc:`How-tos </how_to_list>`