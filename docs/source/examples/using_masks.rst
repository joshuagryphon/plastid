Excluding (masking) regions of the genome
=========================================

It is often desirable to exclude missing or unreliable data from analysis.
For this purpose, `numpy`_ and `SciPy`_ offer masking tools for numerical data
in the modules :mod:`numpy.ma` and :mod:`scipy.stats.mstats`.

In genomics, specific regions of the genome -- e.g. paralogs, transposons and
repetitive sequences --
are prone to yield missing or unreliable data, because reads from
:term:`high-throughput sequencing` experiments will :term:`multimap <multimapping>`
equally well to each instance of the repeated sequence, creating
ambiguity in the data.

In this tutorial we discuss how to exclude or *mask* regions of the genome
from analysis using the :doc:`/test_dataset`. Masking is discussed in the
following sections:

.. contents::
   :local:
   



.. _masking-manual-mask:

Manual masking
--------------

Creating & adding masks
.......................

Like all :term:`genomic features <feature>` in :data:`Plastid`, masks are 
represented as |SegmentChains|, which can be used to mask other |SegmentChain|
or |Transcript| objects.

First, let's load a transcripts, get its CDS, and create masks that cover
the first and last 5 codons:

.. code-block:: python

   >>> import numpy
   >>> from plastid import *
   
   # load transcripts and count data
   >>> alignments = BAMGenomeArray("SRR609197_riboprofile_5hr_rep1.bam",mapping=FivePrimeMapFactory(offset=14))
   >>> transcripts = list(BED_Reader("merlin_orfs.bed",return_type=Transcript))

   # Get coding region using get_cds()
   >>> demo_cds = transcripts[39].get_cds()

   # Create SegmentChains covering the first and last 5 codons. 
   # We'll use these as masks
   >>> start_codon_mask = list(demo_cds.get_subchain(0,15))
   >>> stop_codon_mask  = list(demo_cds.get_subchain(demo_cds.length-15,demo_cds.length))


The |SegmentChains| we just created can be applied as masks using 
:meth:`~plastid.genomics.roitools.SegmentChain.add_masks`:

.. code-block:: python
   
   # Apply the masks to the demo_cds
   >>> demo_cds.add_masks(*start_codon_mask)
   >>> demo_cds.add_masks(*stop_codon_mask)

   # save masks to a BED file
   >>> fout = open("merlin_start_codon_masks.bed","w")
   >>> for mask in start_codon_mask:
   >>>     fout.write(SegmentChain(mask).as_bed())
   >>>
   >>> fout.close()


Uses of masks
.............

After masks are added, we can get a masked count vector by calling
:meth:`~plastid.genomics.roitools.SegmentChain.get_masked_counts`. This method
returns a :class:`numpy.ma.MaskedArray`, rather than a :class:`numpy.ndarray`.
:class:`~numpy.ma.MaskedArray` objects because they contain all the values,
but ignore masked values when performing operations:

.. code-block:: python

   # count reads, excluding those mapping to masked positions
   >>> demo_cds.get_masked_counts(alignments).sum()
   53.0

Calling :meth:`~plastid.genomics.roitools.SegmentChain.get_counts` after adding
masks will still return an *unmasked* :class:`numpy.ndarray`:

.. code-block:: python

   # count all reads
   >>> demo_cds.get_counts(alignments).sum()
   67.0

Masked positions are also excluded from length measurements, if and only if
:meth:`~plastid.genomics.roitools.SegmentChain.get_masked_length` is called:

.. code-block:: python

   >>> demo_cds.masked_length # length, excluding masked nucleotides
   213

   >>> demo_cds.length # total length
   243


We can also retrieve masks that have been added to a |SegmentChain|, either
as a list of |GenomicSegments| or as a |SegmentChain|:

.. code-block:: python

   >>> demo_cds.mask_segments
   [<GenomicSegment merlin:14615-14630 strand='+'>,
    <GenomicSegment merlin:14843-14858 strand='+'>]

   >>> demo_cds.get_masks_as_segmentchain()
   <SegmentChain segments=2 bounds=merlin:14615-14858(+) name=merlin:14615-14630^14843-14858(+)>


.. _masking-mask-files:

Using mask files
----------------
:term:`Mask files <mask file>` annotate genomic regions that should be masked
from analysis. As any annotation file, these can be in any of many formats
(e.g. `BED`_, `BigBed`_, `GFF3`_, or others).


.. _masking-mask-file-interactive

In interactive Python sessions
..............................

:term:`Mask files <mask file>` can be loaded into a |GenomeHash|, a
dictionary-like object that indexes features by their locations in the genome.
To create a |GenomeHash|:

.. code-block:: python

   # get list of masks
   >>> mask_features = list(BED_Reader("merlin_start_codon_masks.bed"))

   # use GenomeHash to index masks
   >>> mask_hash = GenomeHash(mask_features)

We'll retrieve all the masks in `mask_hash` that overlap `demo_cds` by using
it as a dictionary key:

.. code-block:: python

   # find masks
   >>> demo_masks = mask_hash[demo_cds]
   >>> demo_masks
   [<SegmentChain segments=1 bounds=merlin:14615-14630(+) name=merlin:14615-14630(+)>]

   # add to demo_cds
   >>> for mask_chain in demo_masks:
   >>>    demo_cds.add_masks(*mask_chain)

If the :term:`mask file` is very large, it should be converted to an
:term:`indexed file format` such as `BigBed`_ to save memory.

Indexed annotation files can instead be loaded into |BigBedGenomeHash| and
|TabixGenomeHash|, which take advantage of the indexes present in
`BigBed`_ and `tabix`_-compressed files.


.. _masking-mask-file-command-line

In command-line scripts
.......................

:term:`Mask files <mask file>` can be used by :mod:`command-line scripts <plastid.bin>`
using the argument ``--mask_annotation_files``. For example:

.. code-block:: shell

   # create metagene file that excludes regions in mask_file.bed
   $ metagene generate outbase
                       --landmark cds_start \
                       --annotation_files annotation_file.gtf \
                       --mask_annotation_files mask_file.bed \
                       --mask_annotation_format BED


.. _masking-crossmap-script:

Creating a mask file using the |crossmap| script
------------------------------------------------

The |crossmap| script empirically annotates genomic regions that multimap 
under various alignment criteria, and saves these as a  :term:`mask file`.

Algorithm
.........

|crossmap| uses the following approach (adapted from :cite:`Ingolia2009`):

#. A genome is diced into pseudo-reads (:term:`k-mers <k-mer>`) of a given length.
   The length of the pseudo-read is chosen to conservatively approximate the expected
   read length from a :term:`high-throughput sequencing` experiment. So, for a
   :term:`ribosome profiling` experiment that typically produces 27- to 32-mers,
   one might choose `k` to be 25 or 30.

#. The pseudo-reads are realigned to the genome sequence, permitting a user-configurable
   number of mismatches. Again, the number of mismatches should be chosen to conservatively
   reflect the number of mismatches that will be permitted when the sequencing
   data is aligned.

#. The number of times each pseudo-read aligns is counted. When a pseudo-read
   :term:`multimaps <multimapping>` equally well to more than a single location,
   the genomic position that gave rise to that pseudo-read is annotated as
   repetitive under the given value for `k` and number of mismatches.

#. Repetitive regions are saved in `BED`_ format.

Running
.......

Because |crossmap| internally uses `bowtie`_ for alignments, `bowtie`_
must be installed on your system. Once it is, use ``bowtie-build`` to
build an index of your genome. From the terminal:

.. code-block:: shell

   $ bowtie-build merlin_NC006273-2.fa merlin_NC006273-2

   
Then, run the script. We'll use 26-mers and a 12-nucleotide P-site offset,
allowing 2 mismatches during alignment:

.. code-block:: shell

   $ crossmap -k 26 --offset 12 --mismatches 2 \
              merlin_NC006273-2.fa \
              merlin_NC006273-2 \
              merlin_NC006273-2


..


Considerations for large genomes
................................

For large genomes (e.g. vertebrate, plant, or some *very* big amoebas):

 - |crossmap| can require a ton of memory if genome sequence is stored 
   in a fasta file. If |crossmap| maxes out your system's memory, it may
   be terminated by your system before it completes.
   
   Consider converting the file to a `2bit`_ file to save memory and
   avoid this potential problem

 - |crossmap| can take several days to run, especially if mismatches are
   allowed. Consider using ``--mismatches 0`` if you run into this problem 

 - Using more processes (e.g. via ``-p 2``) will speed |crossmap|'s runtime,
   but will increase its memory footprint, as each process will need its own
   memory space to create and align k-mers from chromosomal sequence

 - By default, |crossmap| creates `BED`_ files. Consider converting these to
   `BigBed`_ files will save substantial amounts of time and memory in the future.
   For instructions, see the documentation for `Jim Kent's utilities`_



-------------------------------------------------------------------------------

See also
--------

- :mod:`plastid.genomics.genome_hash`, which includes additional genome hashes
  for various binary or indexed file formats
  
- The |crossmap| script

- :mod:`numpy.ma` and :mod:`scipy.stats.mstats`
  for lists of `numpy`_ and `SciPy`_ functions that operate on 
  :class:`~numpy.ma.MaskedArray` objects
  
- `Jim Kent's utilities`_ for `BigBed`_ conversion.
  
