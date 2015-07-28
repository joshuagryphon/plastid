Change log
==========

All major changes to ``yeti`` are documented here. Version numbers for the project
follow the conventions described in :pep:`440`. After release 1.0, they will 
also follow `Semantic versioning <http://semver.org/>`_.

 .. note::
 
    This project was initially developed internally under the provisional
    name ``genometools``. At that point, we didn't expect to release the project,
    so we didn't checkt os ee whether the naem was free. It turns out, 
    it was not!
    
    When we realized we would release this project extramurally, we provisionaly
    codenamed it ``yeti``, moved some things around the package to more reasonable
    locations, changed some method names and behaviors for consistency, and reset
    the version number. At this same time, we migrated from our old SVN
    repository to a git repository.


yeti [0.1.1] = [2015-07-23]
---------------------------

Added
.....
  - Created & backpopulated changelog
  - Docstrings re-written for user rather than developer focus
  - AssembledFeatureReader
  - Complete first draft of user manual documentation
  - Readthedocs support for documentation

Changed
.......
  - import of scientific packages now simulated using 'mock' during documentation
    builds by Sphinx

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
      - GF3_to_Transcripts
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
  - Merger of make_tophat_juncs, find_juncs, and merge_juncs into one script
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
   
 


        
