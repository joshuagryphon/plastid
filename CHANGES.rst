Change log
==========

All major changes to ``plastid`` are documented here. Version numbers for the project
follow the conventions described in :pep:`440`. After release 1.0, they will 
also follow `Semantic versioning <http://semver.org/>`_. Until then, they roughly
follow `Semantic versioning <http://semver.org/>`_, with a prepended '0.'

 .. note::
 
    This project was initially developed internally under the provisional
    name ``genometools``, and then later under the codename ``yeti``.

Unreleased
----------
Changes here are all backwards-compatible, and are mostly aesthetic changes
aimed at making the library easier to use with less typing

Added
.....
  - Logical imports- the following commonly-used data structures can
    now be directly imported from the parent package ``plastid``, 
    instead of subpackages/submodules:
    
      - ``GenomicSegment``, ``SegmentChain``, and ``Transcript``
      - All GenomeHashes and GenomeArrays
      - All file readers

  - ``VariableFivePrimeMapFactory`` can now be created from static method ``from_file()``,
    so no need to manually parse text files or create dictionaries

  - ``BAMGenomeArray`` can now be initialized with a list of paths to BAM files,
    in addition or instead of a list of ``pysam.AlignmentFiles``

Changed
.......
  - ``add_three_for_stop_codon()`` reimplemented in Cython, resulting in 2-fold speedup.
    Moved from ``plastid.readers.common`` to ``plastid.genomics.roitools`` (though previosu
    import path still works)

Removed
.......
  - Removed deprecated functions ``BED_to_Transcripts()`` and ``BED_to_SegmentChains``,
    for which ``BED_Reader`` serves as a drop-in replacement



plastid [0.3.2] = [2015-10-01]
------------------------------

Changed
.......
  - Important docstring updates: removed outdated warnings and descriptions


plastid [0.3.0] = [2015-10-01]
------------------------------

Changed
.......
  - Cython implementations of ``SegmentChain`` and ``Transcript`` provide
    massive speedups
  - ``Transcript.cds_start``, ``cds_genome_start``, ``cds_end``, ``cds_genome_end``
    are now managed properties and update each other to maintain
    synchrony
  - ``SegmentChain._segments`` and ``SegmentChain._mask_segments`` are now
    read-only

Deprecated
..........
  - Methods ``SegmentChain.get_length()`` and ``SegmentChain.get_masked_length()``
    are replaced by properties ``SegmentChain.length`` and
    ``SegmentChain.masked_length``

Removed
.......
  - ``sort_segments_lexically()`` and ``sort_segmentchains_lexically()``
    removed, because ``GenomicSegment`` and ``SegmentChain`` now sort
    lexically without help


plastid [0.2.3] = [2015-09-23]
------------------------------

Changed
.......
  - Cython implementations of BAM mapping rules now default,
    are 2-10x faster than Python implementations


plastid [0.2.2] = [2015-09-15]
------------------------------
First release under official name!

Added
.....
  - Major algorithmic improvements to internals & command-line scripts

Changed
.......
  - Reimplemented mapping rules and some internals in Cython,
    giving 2-10x speedup for some operations
  - ``GenomicSegment`` now sorts lexically. Properties are read-only


yeti [0.2.1] = [2015-09-06]
---------------------------

Added
.....
  - Support for extended BED formats now in both import & export,
    in command-line scripts and interactively
  - BED Detail format and known ENCODE BED subtypes now automatically parsed
    from track definition lines
  - Created warning classes DataWarning, FileFormatWarning, and ArgumentWarning
  - parallelized `crossmap` script
  - command line support for more sequence file formats; and a sequence argparser

Changed
.......
  - speed & memory optimizations for `cs generate` script, resulting in 90% memory
    reduction on human genome annotation GrCh38.78
  - ditto `metagene generate` script
  - `crossmap` script does not save kmer files unless --save_kmers is given
  - warnings now given at first (instead of every) occurence
  - lazy imports; giving speed improvements to command-line scripts


yeti [0.2.0] = [2015-08-26]
---------------------------
**Big changes,** including some that are **backwards-incompatible.**
We really think these are for the best, because they improve
compatibility with other packages (e.g. pandas) and make
the package more consistent in design & behavior

Added
.....
  - GenomeArray __getitem__ and __setitem__ now can take
    SegmentChains as arguments
  - Mapping functions for bowtie files now issue warnings
    when reads are unmappable
  - support for 2bit files (via twobitreader) and for
    dicts of strings in SegmentChain.get_sequence
  - various warnings added

Changed
.......
  - pandas compatibility: header rows in all output files no longer have starting '#.
    meaning UPDATE YOUR OLD POSITIONS/ROI FILES
  - __getitem__ from GenomeArrays now returns vectors 5' to 3'
    relative to GenomicSegment rather than to genome. This
    is more consistent with user expectations.
  - _get_valid_X methods of SegmentChain changed to _get_masked_X
    for consistency with documentation and with numpy notation

Removed
.......
  - ArrayTable class & tests


yeti [0.1.1] = [2015-07-23]
---------------------------

Added
.....
  - Created & backpopulated changelog
  - Docstrings re-written for user rather than developer focus
  - AssembledFeatureReader
  - Complete first draft of user manual documentation
  - Readthedocs support for documentation
  - GFF3_TranscriptAssembler now also handles features whose subfeatures
    share `ID` attributes instead of `Parent` attributes.

Changed
.......
  - import of scientific packages now simulated using `mock` during documentation
    builds by Sphinx
  - duplicated attributes in GTF2 column 9 are now catenated & returned as a list
    in attr dict. This is outside GTF2 spec, but a behavior used by GENCODE

Fixed
.....
  - Removed bug from :func:`yeti.bin.metagene.do_generate` that extended
    maximal spanning windows past equivalence points in 3' directions.
    Added extra unit test cases to suit it.
  - GenomeHash can again accept GenomicSegments as parameters to __getitem__.
    Added unit tests for this.

Removed
.......
  - Removed deprecated functions, modules, & classes:
      - GenomicFeature
      - BED_to_Transcripts
      - BigBed_to_Transcripts
      - GTF2_to_Transcripts
      - GFF3_to_Transcripts
      - TagAlignReader


yeti [0.1.0] = [2015-06-06]
---------------------------
First internal release of project under new codename, ``yeti``. Reset version 
number.

Added
.....
  - AssembledFeatureReader, GTF2_TranscriptAssembler, GFF3_TranscriptAssembler
  - GTF2/GFF3 token parsers now issue warnings on repeated keys
  - GFF3 token parsers now return 'Parent', 'Alias', 'Dbxref', 'dbxref', and 'Note'
    fields as lists

Changed
.......
  - Package renamed from ``genometools`` to its provisional codename ``yeti``
  - Reset version number to 0.1.0
  - Refactored existing readers to descent from AssembledFeatureReader
  - Migration from old SVN to GIT repo
  - Renaming & moving of functions, classes, & modules for consistency and
    to avoid name clashes with other packages
 
        ==================================  ====================================
        Old name                            New Name
        ----------------------------------  ------------------------------------
        GenomicInterval                     GenomicSegment
        IVCollection                        SegmentChain
        NibbleMapFactory                    CenterMapFactory
        genometools.genomics.ivtools        yeti.genomics.roitools
        genometools.genomics.readers        yeti.readers
        genometools.genomics.scriptlib      yeti.util.scriptlib
        ==================================  ====================================


genometools [0.9.1] - 2015-05-21
--------------------------------

Changed
.......
  - renamed suppress_stdr -> capture_stderr

Added
.....
  - capture_stdout decorator


genometools [0.9.0] - 2015-05-20
--------------------------------

Changed
.......
  - All functions that used GenomicFeatures now use IVCollections instead

Removed
.......
  - GenomicFeature support from GenomeHash subclasses
  - GenomicFeature support from IVCollection and GenomicInterval overlap
    end quality criteria

Deprecated
..........
- GenomicFeature


genometools [0.8.3] - 2015-05-19
--------------------------------

Added
.....
  - Included missing `.positions` and `.sizes` files into egg package

genometools [0.8.2] - 2015-05-19
--------------------------------

Changed
.......
  - Test data now packaged in eggs
  - updated documentation

Fixed
.....
  - Bug in cleanup for test_crossmap
  - Bug in setup.py


genometools [0.8.1] - 2015-05-18
--------------------------------

Added
.....
  - Python 3.0 support
  - Support for tabix-compressed files. Creation of TabixGenomeHash


Changed
.......
  - Propagate various attributes to sub-features (utr_ivc, CDS) from Transcript
  - Propagate all attributes to sub-features during GTF export from Transcript
  - GTF2 export of Transcript objects now generates 'start_codon' and
    'stop_codon' features
  - Update of setup.py and Makefile to make dev vs distribution eggs
  - 'transcript_ids' column of 'cs generate' position file now sorted before
    comma join.


genometools [0.8.2015-05-08] - 2015-05-08
-----------------------------------------

Changed
.......
  - Merger of `make_tophat_juncs`, `find_juncs`, and `merge_juncs` into one script
  - Standardization of column names among various output files:
    region, regions_counted, counts
  - Standardized method names in IVCollection: get_valid_counts, get_valid_length,
    get_length, get_counts, et c
  - IVCollection/Transcript openers/assemblers all return generators and can take
    multiple input files
  - IVCollection/Transcript openers/assemblers return lexically-sorted objects
  - Update to GFF3 escaping conventions rather than URL escaping. Also applied to 
    GTF2 files
  - Refactors to `cs` script, plus garbage collection to reduce memory usage
 
Added
.....
  - Changelog
  - Implementation of test suites
  - Lazy assembly of GFF3 and GTF2 files to save memory in
    `GTF2_TranscriptAssembler` and `GFF3_TranscriptAssembler`
  - BigBed support, creation of BigBedReader and BigBedGenomeHash. AutoSQL support
  - Supported for truncated BED formats
  - P-site offset script
  - `get_count_vectors` script
  - `counts_in_region` script
  - UniqueFifo class
  - Decorators: `parallelize, suppress_stderr, in_separate_process`
  - variableStep export for `BAMGenomeArray`
  - Support of GTF2 "frame" attribute for CDS features

Fixed
.....
  - Bugfixes in minus strand offsets in crossmaps
  - Fixed bug where minus strand crossmap features were ignored
  - Bugfixes in CDS end export from Transcript when CDSes ended at the endpoint
    of internal but not terminal introns on plus-strand transcripts


Deprecated
..........
  - spliced_count_files
  - Ingolia file tagalign import
  - Deprecation of `GTF2_to_Transcripts` and `GFF3_to_Transcripts`
   
 


        
