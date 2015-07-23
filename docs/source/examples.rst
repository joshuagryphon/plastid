Examples
========

The following examples walk through analyses with :data:`yeti` and 
illustrate how to use the libraries. We suggest downloading the
:doc:`test_dataset` and following along.


Cookbook
--------
The following tutorials describe how to perform simple analyses, with minimal
discussion:

 .. toctree::
    :hidden:
    :glob:
    
    /examples/*
    
====================================    ======================================================================================
**Tutorial**                            **Description**
------------------------------------    --------------------------------------------------------------------------------------

:doc:`/examples/gene_expression`        Calculate read densities for :term:`RNA-seq`: & :term:`ribosome profiling`,
                                        and use these to estimate translation efficiency. Then, test for 
                                        differential gene expression.

:doc:`/examples/genome_setup`           Set up a genome for use with :data:`yeti`. Download sequence and
                                        :term:`annotation` files, prepare a :term:`mask file` via |crossmap|, and
                                        discuss other considerations relevant to genomic analyses.

:doc:`/examples/metagene`               Perform a positional :term:`metagene analysis <metagene>`, using :term:`ribosome profiling`
                                        data at the start codon as an example. Then, develop metagene analysis around
                                        a custom landmark for use with other data types.

:doc:`/examples/p_site`                 Determine a :term:`P-site offset` from :term:`ribosome profiling` data
====================================    ======================================================================================

 .. _concepts-index:
 
In-depth
--------
The documents below discuss various issues and concepts in genomics and
:term:`high-throughput sequencing` in depth.


  - :doc:`concepts/coordinates`
  - :doc:`concepts/multimappers`
  - :doc:`concepts/mapping_rules`



 .. toctree::
    :hidden:
    :glob:
    
    /concepts/*
        

See also
--------
:doc:`/tour`
    A brief overview of the important data structures

:doc:`/concepts`
    Concepts & conventions used in genomics



