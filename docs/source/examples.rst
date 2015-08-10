Tutorials
=========
Tutorials are divided into two sections:

  - The :ref:`examples-cookbook` contains step-by-step instructions for
    performing common tasks, with detailed explanations of each step
  
  - :ref:`examples-concepts` includes longer discussions of issues that arise
    in :term:`high-throughput sequencing` and genomics, and includes code
    examples when appropriate. 
  
We suggest downloading the :doc:`test_dataset` and following along.

 
 .. _examples-cookbook:

Cookbook
--------

 .. toctree::
    :hidden:
    :maxdepth: 4
    :glob:
    
    /examples/*
    
    
    
============================================================================    ===========================================================================================
**Tutorial**                                                                    **Contents**
----------------------------------------------------------------------------    -------------------------------------------------------------------------------------------
:doc:`/examples/count_vector`                                                   Retrieve a vector of :term:`high-throughput sequencing` :term:`counts`
                                                                                at each position in a transcript

:doc:`/examples/using_masks`                                                    Exclude specific regions -- for example, repetitive genome sequence that gives rise to
                                                                                :term:`multimapping` reads -- from analysis. Discussion of :term:`mask files <mask file>`.

:doc:`Gene expression & translation efficiency </examples/gene_expression>`     Compute gene expression measurements and translation efficicency using :term:`RNA-seq` & :term:`ribosome profiling`, and prepare
                                                                                data for differential expression analysis

:doc:`Custom genome annotations </examples/make_annotation>`                    Make a custom `BED`_, `BigBed`_, `GTF2`_, or `GFF3`_ file containing custom :term:`features <feature>`.

:doc:`Metagene analysis </examples/metagene>`                                   Perform :term:`metagene analysis <metagene>`, using :term:`ribosome profiling`
                                                                                data at the start codon as an example. Then, develop metagene analysis around
                                                                                a custom landmark for use with other data types

:doc:`Ribosomal P-site offsets </examples/p_site>`                              Determine a :term:`P-site offset` from :term:`ribosome profiling` data

:doc:`/examples/phasing`                                                        Estimate :term:`read phasing (triplet periodicity) <sub-codon phasing>` of :term:`ribosome profiling`
                                                                                data
============================================================================    ===========================================================================================


 .. _examples-concepts:
 
In-depth
--------

===============================================================    ===========================================================================================
**Tutorial**                                                       **Contents**
---------------------------------------------------------------    -------------------------------------------------------------------------------------------
:doc:`/concepts/data`                                              Introduction & discussion to the types of data used in genomics, 
                                                                   and the advantages and disadvantages of their various file formats

:doc:`/concepts/coordinates`                                       Primer on the various coordinate systems used in genomics

:doc:`/concepts/multimappers`                                      Issues arising when and strategies for handling :term:`multimapping` reads

:doc:`/concepts/mapping_rules`                                     In-depth discussion of :term:`mapping rules <mapping rule>`, with code examples
                                                                   of how to write your own :term:`mapping rule` for your own sequencing data type
===============================================================    ===========================================================================================


 .. toctree::
    :hidden:
    :maxdepth: 4
    :glob:
    
    /concepts/*
        
