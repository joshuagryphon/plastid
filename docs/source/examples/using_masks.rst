Excluding (masking) regions of the genome
=========================================

Often, it is important to exclude missing or unreliable data from analysis.
`numpy`_ and `SciPy`_ offer tools to do this with numerical data in the modules
:mod:`numpy.ma` and :mod:`scipy.stats.mstats`.

In genomics, specific regions of the genome -- repetitive sequences --
are prone to yield missing or unreliable data, because reads from
:term:`high-throughput sequencing` experiments will :term:`multimap <multimapping>`
equally well to each repetition of the sequence, creating ambiguity
in the data.

In this tutorial we discuss how to exclude or *mask* regions of the genome
from analysis. We use repetitive sequence as an example, but, any region
could be masked. Masking is discussed in the following sections:

  - :ref:`masking-manual-mask` illustrates the effects of masking
    on tabulation of :term:`feature` :term:`counts` and length
  
  - :ref:`masking-mask-files` describes :term:`mask files <mask file>`,
    and how to use them
    :ref:`in interactive Python sessions <masking-mask-file-interactive>`
    and
    :ref:`in command-line scripts <masking-mask-file-command-line>`.

  - :ref:`masking-crossmap-script` explains how to use the |crossmap|
    script to empirically annotate repetitive regions of the genome

In this tutorial, we use the :doc:`/test_dataset`. You are
encouraged to folow along.


 .. _masking-manual-mask:

Manual masking of regions
-------------------------

Genomic features in :data:`yeti` are represented by |SegmentChains|
and |Transcripts|, which are constructed from zero or more |GenomicSegments|.
Masks, which are features, can also be represented as |SegmentChains|
and |GenomicSegments|. Portions of |SegmentChains| and |Transcripts|
can be masked using their :meth:`~yeti.genomics.roitools.SegmentChain.add_masks`
methods::

    >>> import pysam
    >>> import numpy
    >>> from yeti.genomics.roitools import GenomicSegment
    >>> from yeti.readers.bed import BED_Reader
    >>> from yeti.genomics.genome_array import BAMGenomeArray, VariableFiveprimeMapFactory

    >>> # load transcripts and count data
    >>> offset_dict = 
    >>> alignments = BAMGenomeArray([pysam.Samfile(,"rb")],VariableFivePrimeMapFactory(offset_dict))
    >>> transcripts = list(BED_Reader(open("")))

    >>> #this is ribosome profiling data, so we'll look at a coding region
    >>> demo_cds = transcripts[].get_cds()

    >>> # Now, add masks. We'll mask out the first and last 5 codons
    >>> start_codon_masks = list(demo_cds.get_subchain(0,15))
    >>> stop_codon_masks  = list(demo_cds.get_subchain(demo_cds_length-15,demo_cds_length))
    >>> demo_cds.add_masks(*start_codon_masks)
    >>> demo_cds.add_masks(*stop_codon_masks)

After masks are added, we can get a masked count vector by calling
:meth:`~yeti.genomics.roitools.SegmentChain.get_masked_counts`. This method
returns a :class:`numpy.ma.MaskedArray`, rather than a :class:`numpy.ndarray`.
:class:`~numpy.ma.MaskedArray` objects are useful, because they contain 
all values, but ignore masked values when performing operations::

    >>> # Now, get masked counts and masked length
    >>> demo_cds.get_masked_counts(alignments)

    >>> demo_cds.get_masked_counts(alignments).sum()

Calling :meth:`~yeti.genomics.roitools.SegmentChain.get_masked_counts` after adding
masks will still return an *unmasked* :class:`numpy.ndarray`::

    >>> demo_cds.get_counts(alignments)

    >>> demo_cds.get_counts(alignments).sum()

Masked positions are also excluded from length measurements, if and only if
:meth:`~yeti.genomics.roitools.SegmentChain.get_masked_length` is called::

    >>> demo_cds.get_masked_length() # masked length

    >>> demo_cds.get_length() # unmasked length


We can also retrieve masks that have been added to a |SegmentChain|, either
as a list of |GenomicSegments| or as a |SegmentChain|::

    >>> demo_cds.get_masks()

    >>> demo_cds.get_masks_as_segmentchain()



.

 .. _masking-mask-files:

:term:`Mask files <mask file>`
------------------------------
Mask files are :term:`annotation files <annotation>` that define parts of a
genome, transcriptome, or contig to exclude from analysis. They can be
in any annotation format (e.g. `BED`_, `BigBed`_, `GFF3`_, or others),
and can be used to mask any region, for any reason.


 .. _masking-mask-file-interactive

:mod:`GenomeHashes <yeti.genomics.genome_hash>` and :term:`mask files <mask file>` in interactive Python sessions
.................................................................................................................

:term:`Mask files <mask file>` can be loaded into a |GenomeHash|, which
indexes mask by location in the genome. To create a |GenomeHash|::

    >>> from yeti.genomics.genome_hash import GenomeHash

    >>> # load masks
    >>> mask_features = list(BED_Reader(open()))
    >>> mask_hash = GenomeHash(mask_features)

Then, we can search the |GenomeHash| for relevant masks to apply to features::

    >>> demo_masks = mask_hash[demo_cds]
    >>> demo_masks

    >>> demo_cds.add_masks(*demo_masks)

If the :term:`mask file` is very large, it should be converted to an
:ref:`indexed file format` such as `BigBed`_, or a `tabix`_-compressed file
so that mask features don't need to be held in memory. These formats can be
loaded into |BigBedGenomeHash| and |TabixGenomeHash|, respectively.


 .. _masking-mask-file-command-line

Using :term:`mask files <mask file>` in :mod:`command-line scripts <yeti.bin>`
..............................................................................

:term:`Mask files <mask file>` can be used by :mod:`command-line scripts <yeti.bin>`
if a user supplies the ``--mask_annotation_files`` argument. For example, to 
mask regions when creating a :term:`metagene` window file:

 .. code-block:: shell

    $ metagene generate outbase --landmark cds_start --annotation_files anno_file --mask_annotation_files mask_file --mask_annotation_format BED


 .. _masking-crossmap-script:

Creating a :term:`mask file` of repetitive genome sequence using the |crossmap| script
--------------------------------------------------------------------------------------

The |crossmap| script creates a :term:`mask file` that empirically annotates repetitive
genome sequence, using the following approach (introduced in :cite:`Ingolia2009`):

 #. A genome is diced into pseudo-reads (:term:`k-mers <k-mer>`) of a given length.
    The length of the pseudo-read is chosen to conservatively approximate the expected
    read length from a :term:`high-throughput sequencing` experiment. So, for a
    :term:`ribosome profiling` experiment that typically produces 27- to 32-mers,
    one might choose `k` to be 25 or 30.

 #. The :term:`k-mers <k-mer>` are realigned to the genome sequence, permitting a user-configurable
    number of mismatches. Again, the number of mismatches should be chosen to conservatively
    reflect the number of mismatches that will be permitted when data from the
    :term:`high-throughput sequencing` experiment is aligned.

 #. The number of times each :term:`k-mer` aligns is counted. When a k-mer
    :term:`multimaps <multimapping>` equally well to multiple genomic coordinates,
    the genomic position that gave rise to that :term:`k-mer` is annotated as
    repetitive under the given value for `k` and number of mismatches.

 #. Repetitive regions are saved to a `BED`_ file.


Because |crossmap| internally uses `bowtie`_ for alignments, `bowtie`_
must be installed on your system. Once it is, use ``bowtie-build`` to
build an index of your genome. From the terminal:

 .. code-block:: shell

    $ bowtie-build merlin_NC006273-2.fa merlin_NC006273-2

    
Then, run the script. We'll use 26-mers and a 12-nucleotide P-site offset,
allowing 2 mismatches during alignment:

 .. code-block:: shell

    $ crossmap -k 26 --offset 12 --mismatches 2 merlin_NC006273-2.fa merlin_NC006273-2 merlin_NC006273-2


In this example, the `BED`_ file that is produced is quite small.
But, if it were larger, converting it to a `BigBed`_ file using Jim
Kent's ``bedToBigBed`` would
result in memory savings. For instructions on that conversion, see
the documentation for `Jim Kent's utilities`_.

 .. note::

    For mammalian genomes, |crossmap| can take several days to run,
    especially if mismatches are allowed.


-------------------------------------------------------------------------------

See also
--------

 - Module documentation for :mod:`yeti.genomics.genome_hash`
 - The |crossmap| script
 - Module documentation for :mod:`numpy.ma` and :mod:`scipy.stats.mstats`
   for lists of `numpy`_ and `SciPy`_ functions that operate on 
   :class:`~numpy.ma.MaskedArray` objects
 - `Jim Kent's utilities`_ for `BigBed`_ conversion.
   
