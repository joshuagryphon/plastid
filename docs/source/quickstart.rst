Before starting
===============
To get started with genomic analysis, a few elements are important:

  - Some :doc:`data </concepts/data>` to analyze.
    Frequently, the data consist of:
  
      - A :term:`genome annotation <annotation>`, typically downloaded
        from a consortium like `UCSC`_, `Ensembl`_, `ENCODE`_, `SGD`_,
        or `FlyBase`_.
      
      - A genome sequence that matches the annotation
      
      - Some :term:`high-throughput sequencing` data, in the form of
        :term:`read alignments`, preferably in `BAM`_ format. 
      
        Because good alignment tools -- `bowtie`_, `Tophat`_, and others -- already exist,
        and because choosing alignment parameters involves careful consideration,
        :data:`yeti` does not perform this step for you. For help on performing
        alignments, and a discussion of the subtleties involved, see the
        documentation for the read alignment program you use. 
      
    For the curious, we have created a small :doc:`demo dataset </test_dataset>`,
    used throughout the examples in this documentation. 
  
  - Familiarity with a handful of concepts and conventions. We encourage those
    new to sequencing analysis to check :ref:`examples-concepts` and browse as
    needed.
  
  - Some scientific software for data analysis. See :doc:`/installation`
    for instructions on how to get :data:`yeti` and its dependencies, and
    :doc:`/related` for resources that might be useful.
  
