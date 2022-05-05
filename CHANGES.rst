Change log
==========

Major changes to ``plastid`` are documented here. Version numbers for the
project follow the conventions described in :pep:`440` and
`Semantic versioning 2.0.0 <http://semver.org/>`_.


Unreleased
------------------------------

Put new changes here.


plastid [0.6.0] = [2022-05-05]
------------------------------

This is a maintenance release meant to leave this package in a reasonable state,
as it is only sporadically maintained.


Added
.....

- Dockerization, to further control test environments, and make this release
  release more future-proof. Test environments now include Python 3.6 and 3.9.


Changed
.......

- Bumped minimum requirements to reasonable 2022 standards.

- Upgraded embedded Kent & HTSlib source code

- Clarified licenses


Removed
.......

- Dropped support for Python versions 2.7–3.5. These *might* still run, but
  are no longer tested.

- Deprecated classes: ``BPlusTree``, ``RTree``

- Deprecated methods: of ``SegmentChain.get_length()`` and
  ``SegmentChain.get_masked_length()``


plastid [0.5.1] = [2020-05-20]
------------------------------

- Updates to package metadata and docs


plastid [0.5.0] = [2020-05-20]
------------------------------

Changes are predominantly for maintenance, bugfixes, and streamlining.


Changed
.......

- Iterators rewritten for compatibility with Python 3.7 and Python 3.8 (per
  instructions in :pep:`479` )

- As a result, **while readers in** :mod:`plastid.readers` **are still
  generators, they are no longer their own iterators**. This is a bit of an
  obscure point, but implications are:

   - they can still be cast to lists::

     >>> my_list = list(GTF2_TranscriptAssembler("/path/to/foo.gtf"))

   - they can still be used in ``for`` loops::

     >>> for my_transcript in GTF2_TranscriptAssembler("/path/to/foo.gtf"):
     >>>     # pass

   - they *cannot* be used like this::

     >>> transcripts = GTF2_TranscriptAssembler("/path/to/foo.gtf")
     >>> a = next(transcripts)
     # a will be None!!

    instead, wrap the reader in `iter` if you want to use them this way::

     >>> transcripts = iter(GTF2_TranscriptAssembler("/path/to/foo.gtf"))
     >>> a = next(transcripts)
     # a is now a Transcript


Added
.....

- BioConda support (special thanks to ``@lparsons``)

- Testing streamlined via ``tox``, and additional test environments added


Fixed
.....

- Moved calls to ``open()`` into ``with`` context managers to avoid creation of
  stale filehandles

- Rewrote ``setup.py`` to remove requirement for pre-installation of
  ``cython``, ``numpy``, and ``pysam``, as this is now handled by :mod:`pip`

- Removed references to deprecated ``pandas.DataFrame.sort()``, enabling
  compatibility with ``pandas`` versions above 0.2.0

- Improved compatibility with ``numpy`` versions above 1.31.1

- ``VariableFivePrimeMapFactory.from_file()`` and
  ``StratifiedVariableFivePrimeMapFactory.from_file()`` now work on filenames
  as well as file handles, as they were supposed to

- ``StratifiedVariableFivePrimeMapFactory`` now imported by typing
  ``from plastid import *``

- And others as well


Removed
.......

- ``BigBedReader.custom_fields`` was removed in favor of its non-deprecated
  alias, ``BigBedReader.extension_fields``



plastid [0.4.8] = [2017-04-09]
------------------------------

- Fixed a change in `setup.py` that caused Plastid compilation to fail in
  Macintosh environments. Sorry Mac users!



plastid [0.4.7] = [2017-03-06]
------------------------------

This update is minor compared to the release 0.4.6, and was mainly motivated by
updates, bugfixes, and changes required for compatibility with new versions of
``Pysam``


Added
.....

- Support for ``Pysam`` >= 0.10.0

- ``write_pl_table()`` added as a convenience function

- ``--use_mean`` flag added to ``metagene``

- Warnings / better help text


Fixed
.....

- rounding error in ``get_str_from_rgb()``

- ``PSL_Reader()`` now capable of parsing strands from translated `blat` output

- Fixed bug in header parsing in ``PSL_reader``



plastid [0.4.6] = [2016-05-20]
------------------------------

Highlights

- Support for `BigWig`_ files
- Reimplementation of `BigBed`_ file support
- Simplification of syntax / removal of annoyances in both command-line
  scripts and in infrastructure


Added/Changed
.............

File formats
""""""""""""

- Support for `BigWig`_ files. ``BigWigReader`` reads `BigWig`_ files, and 
  ``BigWigGenomeArray``  handles them conveniently.

- ``BigBedReader`` has been reimplemented using Jim Kent's C library, making
  it far faster and more memory efficient.

- ``BigBedReader.search()`` created to search indexed fields included in BigBed
  files, e.g. to find transcripts with a given `gene_id` (if `gene_id` is included
  as an extension column and indexed). To see which fields are searchable,
  use ``BigBedReader.indexed_fields``


Infrastructure
""""""""""""""

- Simplified file opening. All readers can now take filenames in addition
  to open filehandles. No need to wrap filenames in lists any more.
  For example:
   
   .. code-block:: python

    # old way to open GTF2 file
    >>> data = GTF2_TranscriptAssembler(open("some_file.gtf"))

    # new way. Also works with BED_Reader, GTF2_Reader, GFF3_TranscriptAssembler, and others
    >>> data = GTF2_TranscriptAssembler("some_file.gtf")

    # old way to get read alignments from BAM files
    >>> alignments = BAMGenomeArray(["some_file.bam","some_other_file.bam"])

    # new way
    >>> alignemnts = BAMGenomeArray("some_file.bam","some_other_file.bam")

    # old way to open a tabix-indexed file
    >>> data = BED_Reader(pysam.tabix_iterator(open("some_file.bed.gz"),pysam.asTuple()),tabix=True)

    # new way
    >>> data = BED_Reader("some_file.bed.gz",tabix=True)


  To maintain backward compatibility, the old syntax still works

- ``BAMGenomeArray`` can now use mapping functions that return multidimensional
  arrays. As an example we added ``StratifiedVariableFivePrimeMapFactory``,
  which produces a 2D array of counts at each position in a region (columns),
  stratified by read length (rows).
 
- Reformatted & colorized warning output to improve legibility

- ``read_pl_table()`` convenience function for reading tables written
  by command-line scripts into DataFrames, preserving headers, formatting,
  et c


Command-line scripts
""""""""""""""""""""

- All script output metadata now includes command as executed, for easier
  re-running and record keeping

- Scripts using count files get ``--sum`` flag, enabling users to 
  set effective sum of counts/reads used in normalization and RPKM
  calculations

- ``psite``

   - ``--constrain`` option added to ``psite`` to improve performance on
     noisy or low count data.

   - No longer saves intermediate count files. ``--keep`` option added
     to take care of this.

- ``metagene``

   - Fixed/improved color scaling in heatmap output. Color values are now
     capped at the 95th percentile of nonzero values, improving contrast

   - Added warnings for files that appear not to contain UTRs

   - Like ``psite``, no longer saves intermediate count files. ``--keep``
     option added to take care of this.

- ``phase_by_size`` can now optionally use an ROI file from the 
  ``metagene generate`` subprogram. This improves accuracy in higher
  eukaryotes by preventing double-counting of codons when more than
  one transcript is annotated per gene.

- ``cs chart`` file containing list of genes is now optional. If not given,
  all genes are included in comparisons

- ``reformat_transcripts`` is now able to export extended BED columns 
  (e.g. `gene_id`) if the input data has useful attributes. This particularly
  useful when working with large transcript annotations in GTF2/GFF3 format-
  they can now be exported to BED format, and converted to BigBed foramt,
  allowing random access and low memory usage, while preserving gene-transcript
  relationships.


Fixed
.....

- Version parsing bug in setup script. 

- ``@deprecated`` function decorator now gives ``FutureWarning`` instead
  of ``DeprecationWarning``


Deprecated
..........

- ``--norm_region`` option of ``psite`` and ``metagene`` has been deprecated
  and will be removed in ``plastid`` v0.5. Instead, use ``--normalize_over``,
  which performs the same role, except coordinates are specified relative to the
  landmark of interest, rather than entire window. This change is more
  intuitive to many users, and saves them mental math. If both ``--norm_region``
  and ``--normalize_over`` are specified, ``--normalize_over`` will be used.

- ``BigBedReader.custom_fields`` has been replaced with ``BigBedReader.extension_fields``

- ``BigBedReader.chrom_sizes`` has been replaced with ``BigBedReader.chroms``
  for consistency with other data structures

- ``BPlusTree`` and ``RTree`` classes, which will be removed in ``plastid`` v0.5

 


plastid [0.4.5] = [2016-03-09]
------------------------------

Changes here are mostly under the hood, involving improvements in usability,
speed, stability, compatibility, and error reporting. We also fixed up tools
for developers and added entrypoints for custom mapping rules.


Added
.....

- Users can now control verbosity/frequency of warnings via '-v' or '-q' 
  options! By default there should no long screens of DataWarnings
  when processing Ensembl (or other) GTFs.

- ``--aggregate`` option added to ``psite`` script to improve sensitivity
  for low-count data.

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

- updated plotting tools to fetch color cycles from matplotlib versions >= 1.5
   as well as >= 1.3. This corrected a plotting bug in `cs`.

- :meth:`AnnotationParser.get_genome_hash_from_args` now internally uses 
   GFF3_Reader and GTF2_Reader instead of GFF3_TranscriptAssembler and 
   GTF2_TranscriptAssembler, allowing mask files in GTF2/GFF3 foramts
   to be type-agnostic in command-line scripts

- contig names no longer lost when using 2bit files in `crossmap`

- updates to :mod:`~plastid.bin.psite`
 
   - output header in metagene profiles. Sorry about that 

   - fix compatibility problem with new versions of matplotlib

   - now catches a ``ValueError`` that used to be an ``IndexError``
     in earlier versions of :mod:`numpy`.

- Fixed loss-of-ID bug in :meth:`Transcript.get_cds`


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
  ``plastid.genomics.roitools`` (though previous import path still works)

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

Changed
.......
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

.. note::
 
  This project was initially developed internally under the provisional name
  ``genometools``, and then later under the codename ``yeti``. The current
  name, ``plastid`` will not change. Changelogs from earlier versions 
  appear below.


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
.....a

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
   
 

