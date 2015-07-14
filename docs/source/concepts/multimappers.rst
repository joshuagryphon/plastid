Multimapping
============

Some :term:`read alignments` align equally well to multiple parts of a genome
or transcriptome. This can occur when a read derives from a duplicated gene, 
a transposon, a pseudogene, a segmental duplication, or from a repetitive
sequence like heterochromatin.

In the absence of other information, multimapping to reads cannot be unambiguously
assigned to a single position of origin in a genome or transcriptome. Various
approaches have been developed to handle this, such as:

  - discarding multimappers from alignment, and excluding duplicated genomic
    regions from analysis (as in :cite:`Ingolia2009` and :cite:`Dunn2013`).
 
  - randomly assigning multimappers to all possible places in a genome or
    transcriptome from which they could have arisen (the default behavior
    of the `TopHat`_ aligner)
   
  - using uniquely mapping reads from duplicated genes to determine
    the proportions of multimapping reads that should be assigned
    to each duplicate (TODO: example?)
   
:data:`yeti` is compatible with any of these approaches, but provides
tools specifically for determining which regions of the genome cannot
give rise to uniquely mapping reads (see the |crossmap|) script, and
masking these out of subsequent analyses
(see, for example :meth:`~yeti.genomics.roitools.SegmentChain.add_masks`). 
          

