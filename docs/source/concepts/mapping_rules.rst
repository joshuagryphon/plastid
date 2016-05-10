Read mapping functions
======================

Numerous sequencing assays encode interesting biology in various properties of
read alignments withu nucleotide precision. For example:

 - :term:`ribosome profiling` encodes the positions of ribosomal P-sites as a
   function of both read length and alignment coordinates
   (:cite:`Ingolia2009,Ingolia2011`)

 - Bisulfite sequencing encodes methylation status as C-to-T transitions within
   read alignments

 - :term:`DMS-Seq` (:cite:`Rouskin2014`) and Pseudouridine profiling (:cite:)
   encode unstructured regions of RNA or sites of pseudouridine modification,
   respectively, in the 5' termini fo their read aligments.

:data:`plastid` uses configurable :term:`mapping functions <mapping rule>` to
decode the biology of interest from read alignments, and convert the decoded
data into :class:`arrays <numpy.ndarray>` over regions of interest (such as an
array of the number of ribosomes at each nucleotide in a coding region).

This design enables :term:`plastid`'s tools to operate on sequencing data from
virtually any NGS assay, provided the appropriate mapping function and
parameters. 


Mapping functions are described in the following sections: 

.. contents::
   :local:


.. _mapping-rules-provided:

Mapping functions included in ``plastid``
-----------------------------------------
At present, the following mapping functions are provided, although users are
encouraged to :ref:`write their own mapping functions
<mapping-rules-roll-your-own>` as needed. These include:


*Fiveprime end mapping*
   Each read alignment is mapped to its 5' end, or at a fixed offset (in
   nucleotides) from its 5' end
        
*Variable fiveprime end mapping*
   Each read alignment is mapped at a fixed distance from its 5' end, where the
   distance is determined by the length of the read alignment.
     
   This is used for :term:`ribosome profiling` of yeast (:cite:`Ingolia2009`)
   and mammalian cells (:cite:`Ingolia2011`).

*Stratified variable fiveprime end mapping*
   This multidimensional mapping function behaves as variable fiveprime end
   mapping, and additionally returns arrays that are stratified by read length.
   
   For each region of interest, a 2D array is returned, in which each row
   represents a given read length, and each column a nucleotide position in the
   region of interest. In other words, summing over the columns produces the
   same array that would be given by variable fiveprime end mapping.

*Threeprime end mapping*
   Each read alignment is mapped to its 3' end, or at a fixed offset (in
   nucleotides) from its 3' end.
    
*Entire* or *Center-weighted mapping*
   Zero or more positions are trimmed from each end of the read alignment, and
   the remaining `N` positions in the alignment are incremented by `1/N` read
   counts (so that each read is still counted once, when integrated over its
   mapped length).
     
   This is also used for :term:`ribosome profiling` of *E. coli*
   (:cite:`Oh2011`) and *D. melanogaster* (:cite:`Dunn2013`), and RNA-seq. 

In the image below, the same set of :term:`read alignments` from a
:term:`ribosome profiling` experiment is mapped under various functions. Note
the :term:`start codon peak` and :term:`stop codon peak` that appear when reads
are mapped to specific locations:

.. figure:: /_static/images/mapping_rule_demo.png
   :alt: Ribosome profiling data under different mapping functions
   :figclass: captionfigure
   :width: 1080px
   :height: 683px
    
   **Top**: gene model. **Middle**: alignments of :term:`ribosome footprints`,
   displayed as in the `IGV`_ genome browser without a mapping function.
   **Bottom rows**: :term:`Ribosome footprints` mapped under various mapping
   functions.


.. _mapping-rules-command-line:
 
Choosing mapping functions in command-line scripts
..................................................

Mapping functions may be specified to :mod:`command-line scripts <plastid.bin>`
using the following command-line arguments (note, command-line scripts don't at
present support mapping functions that return multidimensional arrays):

   ======================   ====================================
   **Mapping function**     **Argument**
   ----------------------  ------------------------------------
   Fiveprime                ``--fiveprime``
   
   Fiveprime variable       ``--fiveprime_variable``
   
   Threeprime               ``--threeprime``
   
   Center/entire            ``--center``
   
   ======================  ====================================

The following arguments additionally influence how mapping functions behave:

   ==================== =======================================================
   **Argument**         **Behavior**
   -------------------- -------------------------------------------------------
   ``--offset X``        For ``--fiveprime`` or ``--threeprime``, ``X``
                         is taken to be an integer specifying the offset into
                         the read, at which read alignments should be mapped.
   
                         For ``--fiveprime_variable``, ``X`` is taken to be the
                         filename of a two-column tab-delimited text file, in
                         which first column represents read length or the
                         special keyword `'default'`, and the second column
                         represents the offset from the five prime end at
                         which reads of that length should be mapped.
   -------------------- -------------------------------------------------------
   ``--nibble X``      ``X`` is taken to be the number of bases to trim
                         from each end of the read before mapping.
   ==================== =======================================================

See the documentation for individual :mod:`command-line scripts <plastid.bin>`
for a detailed discussion of their arguments.


.. _mapping-rules-interactive: 
 
Choosing mapping functions in interactive Python sessions
.........................................................

Mapping functions in :data:`plastid` are applied when :term:`read alignments`
are imported. Read alignments are held in data structures called *GenomeArrays*
(see :mod:`plastid.genomics.genome_array`).

Alignments in `BAM`_ format can be imported into a |BAMGenomeArray|. Mapping
functions are set via
:meth:`~plastid.genomics.genome_array.BAMGenomeArray.set_mapping`::

   >>> from plastid.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory, CenterMapFactory

   >>> alignments = BAMGenomeArray(["SRR609197_riboprofile_5hr_rep1.bam"])
   
   >>> # map reads 5 nucleotides downstream from their 5' ends
   >>> alignments.set_mapping(FivePrimeMapFactory(offset=5))

and, the mapping function for a |BAMGenomeArray| can be changed at any time::

   >>> # map reads along entire lengths
   >>> alignments.set_mapping(CenterMapFactory())


Alignments in `bowtie`_ format can be imported into a |GenomeArray|. Because
`bowtie`_ files are not sorted or indexed, mapping functions must be applied
upon import, and cannot be changed afterwards::

   >>> from plastid.genomics.genome_array import GenomeArray, five_prime_map
   
   >>> # map reads 5 nucleotides downstream from their 5' ends
   >>> fiveprime_alignments = GenomeArray()
   >>> fiveprime_alignments.add_from_bowtie(open("some_file.bowtie"),five_prime_map,offset=5)

   >>> # map reads along entire lengths
   >>> entire_alignments = GenomeArray()
   >>> entire_alignments.add_from_bowtie(open("some_file.bowtie"),center_map)


Method names for the various :term:`mapping functions <mapping rule>` appear below:

================================   ==============================================================    ===============================================================================
**Mapping function**               |GenomeArray|, |SparseGenomeArray|                                |BAMGenomeArray|
--------------------------------   --------------------------------------------------------------    -------------------------------------------------------------------------------

Fiveprime                          :func:`~plastid.genomics.genome_array.five_prime_map`             :class:`~plastid.genomics.map_factories.FivePrimeMapFactory`

Fiveprime variable                 :func:`~plastid.genomics.genome_array.variable_five_prime_map`    :class:`~plastid.genomics.map_factories.VariableFivePrimeMapFactory`

Stratified fiveprime variable      not implemented                                                   :class:`~plastid.genomics.map_factories.StratifiedVariableFivePrimeMapFactory`

Threeprime                         :func:`~plastid.genomics.genome_array.three_prime_map`            :class:`~plastid.genomics.map_factories.ThreePrimeMapFactory`

Center/entire                      :func:`~plastid.genomics.genome_array.center_map`                 :class:`~plastid.genomics.map_factories.CenterMapFactory`

================================   ==============================================================    ===============================================================================


.. _mapping-rules-roll-your-own:

Writing new mapping functions
-----------------------------
Mapping functions in :data:`plastid` are implemented as callables (functions or
callable class instances). Mapping functions for |BAMGenomeArray| require the
following signatures:

Parameters
..........
alignments
   A list of :term:`read alignments` represented as
   :class:`pysam.AlignedSegment` objects. These correspond to the alignments
   that will be mapped. Typically, these overlap `segment`.

segment
   A |GenomicSegment| corresponding to a region of interest


Return values
.............
list
   A list of :term:`read alignments` (:class:`pysam.AlignedSegment`) that map
   within `segment` under the mapping function implemented by the function.

:class:`numpy.ndarray`
   An array of values at each position in `segment`, from left-to-right /
   lowest-to-highest coordinates relative to the genome (not relative to the
   segment).
   
   More generally, if the mapping function returns a multi-dimensional array,
   the last dimension must represent the positions in `segment` (e.g., for a 2D
   array, the columns would represent the nucleotide positions).


Example 1: Fiveprime alignment mapping
......................................
This mapping function maps :term:`read alignments` to their 5' ends, allowing an
optional offset::

    >>> import numpy
    >>> import warnings

    >>> def fiveprime_map_function(alignments,segment,offset=0):
    >>>     reads_out = []         
    >>>     count_array = numpy.zeros(len(segment))
    >>>     for read in alignments:
    >>>         if offset > len(read.positions):
    >>>             warnings.warn("Offset %snt greater than read length %snt. Ignoring." % (offset,len(read)),
    >>>                           UserWarning)
    >>>             continue # skip read if offset is outside read boundaries
    >>>             
    >>>         # count offset 5' to 3' if the `segment` is on the plus-strand
    >>>         # or is unstranded
    >>>         if segment.strand == in ("+","."):
    >>>             p_site = read.positions[offset]
    >>>         # count offset from other end if `segment` is on the minus-strand
    >>>         else:
    >>>             p_site = read.positions[-offset - 1]
    >>>          
    >>>         if p_site >= segment.start and p_site < segment.end:
    >>>             reads_out.append(read)
    >>>             count_array[p_site - seg.start] += 1
    >>>             
    >>>    return reads_out, count_array

But, |BAMGenomeArray| will only pass the parameters `alignments` and `segment`
to mapping functions. To specify an offset, use a wrapper function::

   >>> def MyFivePrimeMapFactory(offset=0):
   >>>    def new_func(alignments,segment):
   >>>       return fiveprime_map_function(alignments,segment,offset=offset)
   >>>
   >>>    return new_func

   >>> alignments = BAMGenomeArray(["SRR609197_riboprofile_5hr_rep1.bam"])
   >>> alignments.set_mapping(MyFivePrimeMapFactory(offset=5))   


Example 2: mapping alignments to their mismatches
.................................................
`BAM`_ files contain rich information about read alignments, and these are
exposed to us via :class:`pysam.AlignedSegment`. This mapping function maps
:term:`read alignments` to sites where they mismatch a reference genome,
assuming the alignment contains no indels. Mismatch information is pulled from
the `MD` tag for each read alignment::

   >>> import re
   >>> nucleotides = re.compile(r"[ACTGN]")
   >>> 
   >>> def mismatch_mapping_function(alignments,segment):
   >>>     reads_out = []
   >>>     count_array = numpy.zeros(len(segment))
   >>>     for read in alignments:
   >>>         for tag,val in read.tags:
   >>>             # we are also assuming no indels, which would make parsing MD more complicated.
   >>>             #
   >>>             # mismatches are in stored in `MD` tag of reach alignment in SAM/BAM files
   >>>             # for see MD tag structure http://samtools.sourceforge.net/SAM1.pdf
   >>>             # they basically look like numbers of matches separated by
   >>>             # the letter that mismatches. e.g. 12A15C22
   >>>             # means: 12 matches, followed by mismatch 'A', followed by 15 matches,
   >>>             #        followed by mismatch 'C', followed by 22 matches
   >>>             #
   >>>             # convert MD tag to a vector of positions that mismatch
   >>>             if tag == "MD":
   >>>                 mismatched_positions  = numpy.array([int(X) for X in re.split(nucleotides,val)[:-1]])
   >>>                 mismatched_positions += numpy.arange(len(mismatched_positions))
   >>>     
   >>>         # figure out coordinate of mismatch with respect to genome and `segment`
   >>>         for pos in mismatched_positions:
   >>>             genome_position = read.positions[pos]
   >>>             segment_position = genome_position - segment.start
   >>>             count_array[segment_position] += 1
   >>>     
   >>>     return reads_out, count_array

          
This mapping function may then be used as above::

   >>> alignments.set_mapping(mismatch_mapping_function)      



See also
--------
 - :doc:`/devinfo/entrypoints` for details on making custom mapping functions
   accessible to command-line scripts

 - :doc:`P-site mapping </examples/p_site>` example, in which a mapping function
   for :term:`ribosome profiling` data is derived and applied
    
 - Module documentation for :mod:`plastid.genomics.map_factories` and
   :mod:`plastid.genomics.genome_array`, which provide more details on mapping
   functions, |BAMGenomeArrays|, and |GenomeArrays|

