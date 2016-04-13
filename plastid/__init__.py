#!/usr/bin/env python
"""Welcome to plastid!

This package contains various utilities for analyzing high-throughput sequencing
data, with an emphasis on simplicity for users. To this end, this package provides:

  #. A set of command-line scripts that implement common sequencing workflows
     (see |bin|).
  
  #. Readers that abstract data from various file formats into a minimal set of
     object types. These object types define APIs that easily interface with
     existing scientific tools, such as the `SciPy`_ stack (see |genomics| and
     |readers|)

  #. Tools to facilitate writing command-line scripts (see |scriptlib|)


Package overview
----------------
plastid is divided into the following subpackages:

    ==============    =========================================================
    Package           Contents
    --------------    ---------------------------------------------------------
    |bin|             Command-line scripts
    |genomics|        Classes and functions to manipulate genome annotations, alignments, and quantitative data
    |plotting|        Tools for plotting
    |readers|         Parsers for various file formats
    |util|            Utilities (e.g. function decorators, exceptions, argument parsers)
    |test|            Unit and functional tests (requires download of test datasets)
    ==============    =========================================================
     
"""
__version__ = "0.4.6"
__author__  = "Joshua Griffin Dunn"
import matplotlib
matplotlib.use("agg")

from plastid.genomics.roitools import GenomicSegment, SegmentChain, Transcript
from plastid.genomics.genome_array import BAMGenomeArray, BigWigGenomeArray, GenomeArray, SparseGenomeArray, variable_five_prime_map, five_prime_map, center_map, three_prime_map
from plastid.genomics.genome_hash import GenomeHash, TabixGenomeHash, BigBedGenomeHash
from plastid.genomics.map_factories import VariableFivePrimeMapFactory, FivePrimeMapFactory, CenterMapFactory, ThreePrimeMapFactory, SizeFilterFactory

from plastid.readers.bed import BED_Reader
from plastid.readers.bigbed import BigBedReader
from plastid.readers.gff import GTF2_Reader, GFF3_Reader, GTF2_TranscriptAssembler, GFF3_TranscriptAssembler
from plastid.readers.psl import PSL_Reader
from plastid.readers.bigwig import BigWigReader

from plastid.util.io.openers import read_pl_table

from plastid.util.services.exceptions import formatwarning

