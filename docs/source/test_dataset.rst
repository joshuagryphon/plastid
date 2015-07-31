Demo dataset
============

 .. TODO update the demo dataset filename to include package name

:download:`Download the demo dataset here <_static/demo.tar.bz2>`.

We have put together a small demo dataset that is used in the :doc:`tour`
and in :doc:`examples`. It consists of feature :term:`annotations`, 
:term:`ribosome profiling data`, and :term:`RNA-seq` from the merlin 
(laboratory) strain of human cytomegalovirus (hCMV).

The following files are included:

    ======================================================  =======================================================================  ============================================
    **Filename**                                            **Contents**                                                             **Source**
    ------------------------------------------------------  -----------------------------------------------------------------------  --------------------------------------------
    ``merlin_NC006273-2.fa``                                Sequence of hCMV merlin strain                                           `NCBI`_

    ``merlin_orfs.bed``, ``merlin_orfs.gtf``                Coding region models for hCMV strain, plus estimated UTRs                :cite:`Stern-Ginossar2012` (CDS).
                                                                                                                                     5' UTRs estimated as 50 nt upstream of CDS. 3' UTRs estimated as 100 nt downstream of CDS. 

    ``SRR609197_riboprofile.bam``                           :term:`Ribosome profilingi` data, 5 hours post hCMV infection            :cite:`Stern-Ginossar2012`

    ``SRR592963_rnaseq.bam``                                :term:`RNA-seq` data, 5 hours posth CMV infection                        :cite:`Stern-Ginossar2012`
    ======================================================  =======================================================================  ============================================

