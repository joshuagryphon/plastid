Change log
==========

All major changes to ``plastid`` are documented here. Version numbers for the
project follow the conventions described in :pep:`440`. After release 1.0, they
will also follow `Semantic versioning <http://semver.org/>`_. Until then, they
roughly follow `Semantic versioning <http://semver.org/>`_, with a prepended
'0.'

  .. note::
  
     This project was initially developed internally under the provisional name
     ``genometools``, and then later under the codename ``yeti``. The current
     name, ``plastid`` will not change.



Unreleased
----------

Changes here are mostly under the hood, involving improvements in speed,
stability, compatibility, and error reporting. In addition, entrypoints
were created to allow custom mapping rules to be used from command-line
scripts.


Added
.....

  - Users can now control verbosity/frequency of warnings via '-v' or '-q' 
    options! By default there should no long screens of DataWarnings
    when processing Ensembl (or other) GTFs.

  - Created entrypoints for allowing users to use custom mapping rules
    in the command line scripts:

      - ``plastid.mapping_rules`` for specifying new mapping functions
      - ``plastid.mapping_options`` for specifying any other command-line
        arguments they consume
  
    Detailed instructions for use in the *developer info* section
    of `<plastid.readthedocs.org>`_.

  - Argument parsing classes that replace methods deprecated below:
  
      - :class:`~plastid.util.scriptlib.argparsers.AlignmentParser`
      - :class:`~plastid.util.scriptlib.argparsers.AnnotationParser`
      - :class:`~plastid.util.scriptlib.argparsers.MaskParser`
      - :class:`~plastid.util.scriptlib.argparsers.SequenceParser`
      - :class:`~plastid.util.scriptlib.argparsers.PlottingParser`


Fixed
.....

  - :mod:`~plastid.bin.psite` now catches a ``ValueError`` that used to be
    an ``IndexError`` in earlier versions of :mod:`numpy`.

  - updated plotting tools to fetch color cycles from matplotlib versions >= 1.5
    as well as >= 1.3. This corrected a plotting bug in `cs`.

  - :meth:`AnnotationParser.get_genome_hash_from_args` now internally uses 
    GFF3_Reader and GTF2_Reader instead of GFF3_TranscriptAssembler and 
    GTF2_TranscriptAssembler, allowing GTF2/GFF3 masks to be type-agnostic

  - bug lossing in result of contig names when using 2bit files in `crossmap`

  - updates to `psite`
  
      - output header in metagene_profiles. Sorry about that 

      - fix compatibility problem with new versions of matplotlib


Changed
.......

  - :func:`~plastid.util.services.decorators.deprecated` function decorator
    now optionally takes parameters indicating the future version of plastid
    in which deprecated features will be removed, and what replacement to use
    instead


Deprecated
..........

  - Argument parsing methods:
  
      - ``get_alignment_file_parser()`` & ``get_genome_array_from_args()``.
        Use :class:`~plastid.util.scriptlib.argparsers.AlignmentParser` instead.
      - ``get_annotation_file_parser()`` & ``get_transcripts_from_args()``,
        ``get_segmentchain_file_parser()`` & ``get_segmentchains_from_args()``
        Use :class:`~plastid.util.scriptlib.argparsers.AnnotationParser` instead.
      - ``get_mask_file_parser()`` & ``get_genome_hash_from_mask_args()``.
        Use :class:`~plastid.util.scriptlib.argparsers.MaskParser` instead.
      - ``get_sequence_file_parser()`` & ``get_seqdict_from_args()``.
        Use :class:`~plastid.util.scriptlib.argparsers.SequenceParser` instead
      - ``get_plotting_parser()``, ``get_figure-from_args()``, & ``get_colors_from_args``.
        Use :class:`~plastid.util.scriptlib.argparsers.PlottingParser` instead
      



plastid [0.4.4] = [2105-11-16]
------------------------------

Although the list of changes is short, this release includes dramatic reductions
in memory usage and speed improvements, as well as a few bug fixes. We recommend
everybody upgrade

Added
.....
  - Fast ``merge_segments()`` function in ``roitools`` module.


Changed
.......
  - 10-100 fold reduction in memory consumed by ``SegmentChain`` objects,
    ``GTF2_TranscriptAssembler`` and ``GFF3_TranscriptAssembler``.  All
    position & mask hashes now lazily evaluated
  - 50-fold fold Speed boosts in ``SegmentChain.overlaps()``,
    ``SegmentChain.covers()`` and other methods for comparing ``SegmentChain``
    and ``Transcript`` objects
  - ``GenomicSegment`` is now hashable, e.g. can be used in sets or dict keys 

Fixed
.....
  - Track naming bug in ``make_wiggle``
  - init bug in ``GenomeHash``



plastid [0.4.3] = [2015-10-28]
------------------------------

Fixed
.....
  - Fixed bug in ``crossmap`` script when run on 2bit files



plastid [0.4.2] = [2015-10-22]
------------------------------

No change in codebase vs 0.4.0. Updated required matplotlib version to 1.4.0. 
Made some changes in sphinx doc config for readthedocs.org, which is still
at matplotlib 1.3.0.



plastid [0.4.0] = [2015-10-21]
------------------------------

This release primarily focuses on ease of use: mainly, it is a lot easier
to do things with fewer lines of code. Imports have been shortened, plotting
tools have been added, and scripts now produce more informative output.


Added
.....
   - Logical imports: the following commonly-used data structures can now be
     directly imported from the parent package ``plastid``, instead of
     subpackages/submodules:
     
       - ``GenomicSegment``, ``SegmentChain``, and ``Transcript``
       - All GenomeHashes and GenomeArrays
       - All file readers

   - ``VariableFivePrimeMapFactory`` can now be created from static method
     ``from_file()``, so no need to manually parse text files or create
     dictionaries

   - ``BAMGenomeArray`` can now be initialized with a list of paths to BAM
     files, in addition or instead of a list of ``pysam.AlignmentFiles``

   - **Plotting improvements**

       - ``plastid.plotting`` package, which includes tools for making MA plots,
         scatter plots with marginal histograms, metagene profiles, et c

       - more informative plots made in ``metagene``, ``psite``,
         ``phase_by_size``, and ``cs`` scripts

       - support for matplotlib stylesheets, colormaps, et c in all command-line
         scripts


Changed
.......
   - ``add_three_for_stop_codon()`` reimplemented in Cython, resulting in 2-fold
     speedup.  Moved from ``plastid.readers.common`` to
     ``plastid.genomics.roitools`` (though previosu import path still works)

Fixed
.....
   - Fixed IndexError in ``psite`` that arose when running with the latest
     release of numpy, when generating a read profile over an empty array

   - Legends/text no longer get cut off in plots

Removed
.......
   - Removed deprecated functions ``BED_to_Transcripts()`` and
     ``BED_to_SegmentChains``, for which ``BED_Reader`` serves as a drop-in
     replacement



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
   - ``Transcript.cds_start``, ``cds_genome_start``, ``cds_end``,
     ``cds_genome_end`` are now managed properties and update each other to
     maintain synchrony
   - ``SegmentChain._segments`` and ``SegmentChain._mask_segments`` are now
     read-only

Deprecated
..........
   - Methods ``SegmentChain.get_length()`` and
     ``SegmentChain.get_masked_length()`` are replaced by properties
     ``SegmentChain.length`` and ``SegmentChain.masked_length``

Removed
.......
   - ``sort_segments_lexically()`` and ``sort_segmentchains_lexically()``
     removed, because ``GenomicSegment`` and ``SegmentChain`` now sort lexically
     without help


plastid [0.2.3] = [2015-09-23]
------------------------------

Changed .......
   - Cython implementations of BAM mapping rules now default, are 2-10x faster
     than Python implementations


plastid [0.2.2] = [2015-09-15]
------------------------------

First release under official name!

Added
.....
   - Major algorithmic improvements to internals & command-line scripts

Changed
.......
   - Reimplemented mapping rules and some internals in Cython, giving 2-10x
     speedup for some operations
   - ``GenomicSegment`` now sorts lexically. Properties are read-only


yeti [0.2.1] = [2015-09-06]
---------------------------

Added
.....
   - Support for extended BED formats now in both import & export, in
     command-line scripts and interactively
   - BED Detail format and known ENCODE BED subtypes now automatically parsed
     from track definition lines
   - Created warning classes DataWarning, FileFormatWarning, and ArgumentWarning
   - parallelized `crossmap` script
   - command line support for more sequence file formats; and a sequence
     argparser

Changed
.......
   - speed & memory optimizations for `cs generate` script, resulting in 90%
     memory reduction on human genome annotation GrCh38.78
   - ditto `metagene generate` script
   - `crossmap` script does not save kmer files unless --save_kmers is given
   - warnings now given at first (instead of every) occurence
   - lazy imports; giving speed improvements to command-line scripts


yeti [0.2.0] = [2015-08-26]
---------------------------

**Big changes,** including some that are **backwards-incompatible.** We
really think these are for the best, because they improve compatibility
with other packages (e.g. pandas) and make the package more consistent
in design & behavior

Added
.....
   - GenomeArray __getitem__ and __setitem__ now can take SegmentChains as
     arguments
   - Mapping functions for bowtie files now issue warnings when reads are
     unmappable
   - support for 2bit files (via twobitreader) and for dicts of strings in
     SegmentChain.get_sequence
   - various warnings added

Changed
.......
   - pandas compatibility: header rows in all output files no longer have
     starting '#.  meaning UPDATE YOUR OLD POSITIONS/ROI FILES
   - __getitem__ from GenomeArrays now returns vectors 5' to 3' relative to
     GenomicSegment rather than to genome. This is more consistent with user
     expectations.
   - _get_valid_X methods of SegmentChain changed to _get_masked_X for
     consistency with documentation and with numpy notation

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
   - GFF3_TranscriptAssembler now also handles features whose subfeatures share
     `ID` attributes instead of `Parent` attributes.

Changed
.......
   - import of scientific packages now simulated using `mock` during
     documentation builds by Sphinx
   - duplicated attributes in GTF2 column 9 are now catenated & returned as a
     list in attr dict. This is outside GTF2 spec, but a behavior used by
     GENCODE

Fixed
.....
   - Removed bug from :func:`yeti.bin.metagene.do_generate` that extended
     maximal spanning windows past equivalence points in 3' directions.  Added
     extra unit test cases to suit it.
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
   - GFF3 token parsers now return 'Parent', 'Alias', 'Dbxref', 'dbxref', and
     'Note' fields as lists

Changed
.......
   - Package renamed from ``genometools`` to its provisional codename ``yeti``
   - Reset version number to 0.1.0
   - Refactored existing readers to descent from AssembledFeatureReader
   - Migration from old SVN to GIT repo
   - Renaming & moving of functions, classes, & modules for consistency and to
     avoid name clashes with other packages
  
         ==================================  ====================================
         Old name                            New Name
         ----------------------------------  ------------------------------------
         GenomicInterva                      GenomicSegment
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
   - GenomicFeature support from IVCollection and GenomicInterval overlap end
     quality criteria

Deprecated
..........
- GenomicFeature


genometools [0.8.3] - 2015-05-19
--------------------------------

Added .....
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

Added .....
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
   - Merger of `make_tophat_juncs`, `find_juncs`, and `merge_juncs` into one
     script
   - Standardization of column names among various output files: region,
     regions_counted, counts
   - Standardized method names in IVCollection: get_valid_counts,
     get_valid_length, get_length, get_counts, et c
   - IVCollection/Transcript openers/assemblers all return generators and can
     take multiple input files
   - IVCollection/Transcript openers/assemblers return lexically-sorted objects
   - Update to GFF3 escaping conventions rather than URL escaping. Also applied
     to GTF2 files
   - Refactors to `cs` script, plus garbage collection to reduce memory usage
  
Added
.....
   - Changelog
   - Implementation of test suites
   - Lazy assembly of GFF3 and GTF2 files to save memory in
     `GTF2_TranscriptAssembler` and `GFF3_TranscriptAssembler`
   - BigBed support, creation of BigBedReader and BigBedGenomeHash. AutoSQL
     support
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
    
  


         
