Read phasing in :term:`ribosome profiling`
==========================================

Because the ribosome takes 3-nucleotide steps along the coding region
it is translating, :term:`read alignments` from :term:`ribosome profiling` 
data frequently exhibit :term:`triplet periodicity` when they are mapped
to their :term:`P-site offsets <P-site offset>` (:cite:`Ingolia2009`):

 .. TODO: phasing figure

 .. figure: 
    :alt: Ribosome phasing genome browser examples
    :figclass: captionfigure

    Triplet periodicity across a coding region (source: )

This periodicity is a useful feature, as it can be used to infer the reading
frame in which a coding region is decoded:

 .. TODO: insert phasing chart figure

 .. figure:
    :alt: Phasing differs between reading frames
    :figclass: captionfgure

    Phasing provides unique signatures of reading frames

In this tutorial, we show how to measure :term:`sub-codon phasing` using
|phase_by_size|.

 .. note::

    Heavily biased nucleases -- for example, micrococcal
    nuclease, used for ribosome profiling *E. coli* (:cite:`Oh2011`) and
    *D. melanogaster* (:cite:`Dunn2014`) -- introduce uncertainty into
    P-site mapping, which obscures :term:`phasing <sub-codon phasing>`
    in these datasets.


Measuring :term:`phasing <sub-codon phasing>` using the |psite| script
......................................................................
 .. TODO: fill out p-site usage section


-------------------------------------------------------------------------------

See also
--------
  - :doc:`/examples/p_site` for a discussion on how to determine the
    :term:`P-site offsets <P-site offset>` to use for a given
    :term:`ribosome profiling` dataset.
  - :doc:`/concepts/mapping_rules` for a discussion on how to apply
    :term:`P-site offsets <P-site offset>` or other mapping rules
