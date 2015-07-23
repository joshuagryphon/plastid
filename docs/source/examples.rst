Tutorials
=========
Tutorials are divided into two sections:

  - The :ref:`examples-cookbook` contains step-by-step instructions for
    performing common tasks, with detailed explanations of each operation
  
  - :ref:`examples-concepts` includes detailed discussions of issues that arise
    in :term:`high-throughput sequencing` and genomics, and includes code
    examples when appropriate. 
  
We suggest downloading the :doc:`test_dataset` and following along.

 .. TODO : figure out formatting for this page
 
 .. _examples-cookbook:

Cookbook
--------
Some blurb

 .. toctree::
    :hidden:
    :glob:
    
    /examples/*
    
====================================    ======================================================================================
**Tutorial**                            **Description**
------------------------------------    --------------------------------------------------------------------------------------

:doc:`/examples/gene_expression`        Calculate read densities for :term:`RNA-seq` & :term:`ribosome profiling`,
                                        and use these to estimate translation efficiency. Then, test for 
                                        differential gene expression.

:doc:`/examples/using_masks`            Annotate regions of a genome (or gene or transcript) that are difficult to analyze
                                        due to repetitive or low-complexity sequence, and exclude these from analysis

:doc:`/examples/metagene`               Perform a positional :term:`metagene analysis <metagene>`, using :term:`ribosome profiling`
                                        data at the start codon as an example. Then, develop metagene analysis around
                                        a custom landmark for use with other data types.

:doc:`/examples/p_site`                 Determine a :term:`P-site offset` from :term:`ribosome profiling` data
====================================    ======================================================================================


 .. _examples-concepts:
 
In-depth
--------
Some other blurb
   - :doc:`concepts/coordinates`
   - :doc:`concepts/multimappers`
   - :doc:`concepts/mapping_rules`


 .. toctree::
    :hidden:
    :glob:
    
    /concepts/*
        
