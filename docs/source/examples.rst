Examples
========

The following examples walk through analyses with :data:`yeti` and 
illustrate how to use the libraries. We suggest those new to 
sequencing or bioinformatics download the :doc:`test_dataset` and
follow along.

:doc:`examples/gene_expression`
    In this example, we calculate mRNA expression, ribosome footprint denstiy,
    and translation efficiency in a joint :term:`RNA-seq` / :term:`ribosome profiling`
    experiment. In addition, we test for differential gene expression at transcriptional
    and translational levels.

:doc:`examples/genome_setup`
    In this example, we download sequence and annotation files, prepare a
    :term:`mask file` via the |crossmap| script, and discuss other considerations
    relevant to genomic analyses.

:doc:`examples/metagene`
    In this example, we perform :term:`metagene analysis` of
    :term:`ribosome profiling` data surrounding a start codon. We then describe
    how to perform custom metagene analyses on any data type, and implement a
    functions to perform metagene analyses around the largest peak of 
    ribosome footprint density in a coding region.

:doc:`examples/p_site`
    In this document, we demonstrate how to determine a :term:`P-site offset`
    from :term:`ribosome profiling` data.

 .. toctree::
    :hidden:

    examples/gene_expression
    examples/genome_setup
    examples/metagene
    examples/p_site


See also
--------
:doc:`tour`
    A brief overview of the important data structures

:doc:`concepts`
    Concepts & conventions used in genomics



