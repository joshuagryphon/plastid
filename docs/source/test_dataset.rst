Demo dataset
============

We have put together a small demo dataset that is used in the :doc:`tour`
and in :doc:`examples`. It consists of feature :term:`annotations`, 
:term:`ribosome profiling data`, and :term:`RNA-seq` from the merlin 
(laboratory) strain of human cytomegalovirus (hCMV).


Downloads:

 - `demo dataset part one`_, for all of the tutorials
 
 - `demo dataset part two`_, specifically used in
   :doc:`/examples/gene_expression`.



The following files are included:

   ======================================================  =======================================================================  ====================================================================================================================================
   **Filename**                                            **Contents**                                                             **Source**
   ------------------------------------------------------  -----------------------------------------------------------------------  ------------------------------------------------------------------------------------------------------------------------------------
   ``merlin_NC006273-2.fa``                                Sequence of hCMV merlin strain                                           `GenBank, accession no. NC_006273.2 <http://www.ncbi.nlm.nih.gov/nuccore/NC_006273.2>`_

   ``merlin_orfs.bed``, ``merlin_orfs.gtf``                Coding region models for hCMV strain, plus estimated UTRs                :cite:`Stern-Ginossar2012` (CDS).
                                                                                                                                    5' UTRs estimated as 50 nt upstream of CDS. 3' UTRs estimated as 100 nt downstream of CDS. 

   ``SRR609197_riboprofile_rep1_5hr.bam``                  :term:`Ribosome profiling` data, 5 hours post hCMV infection,            :cite:`Stern-Ginossar2012`,
                                                           aligned to hCMV merlin strain genome sequence                            raw data available at `SRA, accession no. SRR609197 <http://www.ncbi.nlm.nih.gov/sra/?term=SRR609197>`_

   ``SRR592963_rnaseq_rep1_5hr.bam``                       :term:`RNA-seq` data, 5 hours post CMV infection,                        :cite:`Stern-Ginossar2012`,
                                                           aligned to hCMV merlin strain genome sequence                            raw data available at `SRA, accession no. SRR592963 <http://www.ncbi.nlm.nih.gov/sra/?term=SRR592963>`_
   ======================================================  =======================================================================  ====================================================================================================================================


The additional dataset includes further replicates, as well as timepoint data
from 24 hours post-infection.
