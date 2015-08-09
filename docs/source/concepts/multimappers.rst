Ambiguous read alignments
=========================

Some high-throughput sequencing reads align equally well (:term:`multimap <multimapping>`)
to multiple parts of a genome or transcriptome. This can occur when a read derives
from repeated sequence, such as a duplicated gene, transposon, or pseudogene;
or from repetitive sequence like telomeres or heterochromatin.

In the absence of other information, :term:`multimapping` reads cannot
unambiguously be assigned to a single position of origin. Various approaches
have been developed to handle this:

  - discarding :term:`multimappers <multimapping>` from alignment, and
    excluding duplicated genomic regions from analysis using a 
    :term:`mask file` (as in :cite:`Ingolia2009` and :cite:`Dunn2013`)
 
  - counting each copy of repetitive sequence (e.g. each copy of a transposon)
    as a single entity and summing all read alignments across all copies before
    calculating read density

   - randomly assigning each multimapper to one of the possible places in a genome
    or transcriptome from which they could have arisen (the default behavior
    of the `TopHat`_ aligner)
   
  - using uniquely mapping reads surrounding each copy of repeated sequence
    to determine the proportions of multimapping reads that should be assigned
    to each copy (TODO: example?)

  
:data:`yeti` is compatible with any of these approaches, but provides
tools specifically for determining which regions of the genome cannot
give rise to uniquely mapping reads (see the |crossmap|) script, and
for excluding such regions from subsequent analyses
(see, for example :meth:`~yeti.genomics.roitools.SegmentChain.add_masks`). 
          

See also
--------
  - :doc:`/examples/using_masks` for a discussion of masks and
    :term:`mask files <mask file>`
  - the |crossmap| script, which can generate :term:`mask files <mask file>`
    of repetitive sequence, given a genome annotation
