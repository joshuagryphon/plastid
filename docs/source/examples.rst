Examples
========

The following examples walk through analyses with :data:`yeti` and 
illustrate how to use the libraries. We suggest downloading the
:doc:`test_dataset` and following along.

====================================    ======================================================================================
**Example**                             **Description**
------------------------------------    --------------------------------------------------------------------------------------

:doc:`/examples/gene_expression`        Calculate read densities for :term:`RNA-seq`: & :term:`ribosome profiling`,
                                        and use these to estimate translation efficiency. Then, test for 
                                        differential gene expression.

:doc:`/examples/genome_setup`           Set up a genome for use with :data:`yeti`. Download sequence and
                                        :term:`annotation` files, prepare a :term:`mask file` via |crossmap|, and
                                        discuss other considerations relevant to genomic analyses.

:doc:`/examples/metagene`               Perform a positional :term:`metagene analysis`, using :term:`ribosome profiling`
                                        data at the start codon as an example. Then, develop metagene analysis around
                                        a custom landmark for use with other data types.

:doc:`/examples/p_site`                 Determine a :term:`P-site offset` from :term:`ribosome profiling` data
====================================    ======================================================================================


 .. toctree::
    :hidden:

    /examples/gene_expression
    /examples/genome_setup
    /examples/metagene
    /examples/p_site


See also
--------
:doc:`/tour`
    A brief overview of the important data structures

:doc:`/concepts`
    Concepts & conventions used in genomics



