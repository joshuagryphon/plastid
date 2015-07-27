Gene expression analysis
========================

This tutorial illustrates how to measure read density over regions. As 
an example, we look at gene expression (in :term:`raw read counts <counts>` and :term:`RPKM`)
using :term:`RNA-seq` and :term:`ribosome profiling` data. However, the
analysis below could be applied with minimal modification to other 
types of :term:`high-throughput sequencing` data (e.g. ClIP-SEQ, :term:`DMS-Seq`,
et c).

We will do this two ways:

 #. Using the :mod:`~yeti.bin.counts_in_region` script to
    :ref:`count expression automatically <gene-expression-scripts>`.

 #. :ref:`Manually calculating gene expression <gene-expression-interactive>`
    in an interactive Python session

The examples below use the :doc:`/test_dataset` we have assembled. 


Tabulating :term:`read counts <counts>` & :term:`RPKM`
------------------------------------------------------

 .. _gene-expression-scripts:

Via command-line scripts
........................

:data:`yeti` includes two scripts for measuring gene expression:

  * :mod:`~yeti.bin.cs`, which pre-processes a genome anntation and makes
    various heuristic corrections to gene boundaries

  * :mod:`~yeti.bin.counts_in_region`, which does not.

The differences between the scripts are further explained in
:ref:`faq-cs-vs-counts-in-region`.
Here we will use :mod:`~yeti.bin.counts_in_region`.

Our first dataset is :term:`ribosome profiling`, and we will map the ribosomal
P-site at 15 nucleotides from the 3' end of each read (see :cite:`Ingolia2009`).
To specify this, we use the options ``--threeprime --offset 15``.

The data we want to count is in the file ``SRR1562907_chrI.bam``, which we pass
via ``--count_files``. The genes we are interested in counting in this example
are on chromosome I, in the annotation file ``sgd_plus_utrs_chrI.gtf``. Finally,
we will tell the script to save the output in ``riboprofile.txt``.

Putting this together, the script is run from the terminal as:

 .. code-block:: shell

    $ counts_in_region riboprofile.txt --count_files SRR1562907_chrI.bam --annotation_files sgd_plus_utrs_chrI.gtf --threeprime --offset 15

:mod:`~yeti.bin.counts_in_region` will create a tab-delimited text file called
``riboprofile.txt`` containing the results. For detailed documentation of the output
and command-line arguments, see the module documentation for :mod:`~yeti.bin.counts_in_region`.


 .. _gene-expression-interactive:

Manually
........

Even though we have scripts to do this, gene expression can be calculated easily
in an interactive Python session, and it is illustrative to do so. In addition
to caclulating gene expression over entire transcripts, we will also calculate
expression separately over 5' UTRs, coding regions, and 3' UTRs.

First, we need to import a few things::

    >>> import copy
    >>> import pysam
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from yeti.readers.gff import GTF2_TranscriptAssembler
    >>> from yeti.genomics.genome_array import BAMGenomeArray, ThreePrimeMapFactory, CenterMapFactory

.. TODO find an RNA-seq dataset

Then, we'll open our data, storing each dataset in a |BAMGenomeArray|::

    >>> my_datasets = { "ribosome_profiling" : "SRR1562907_chrI.bam",
    >>>                 "RNA-seq"            : "",
    >>>               }

    >>> my_datasets = { K : BAMGenomeArray([pysam.Samfile(V)]) for K,V in my_datasets.items() }

 
Next, we tell the |BAMGenomeArrays| which :term:`mapping rule` to use. We
will map the ribosome-protected footprints to their P-sites, which we estimate
as 15 nucleotides from the 3' end of each read::

    >>> my_datasets["ribosome_profiling"].set_mapping(ThreePrimeMapFactory(offset=15))

We will map the RNA-seq data along the entire length of each read alignment.
Each position in each alignment will be attributed :math:`1.0 / \ell`, where 
:math:`\ell` is the length of the read alignment.
:func:`~yeti.genomics.genome_array.CenterMapFactory` can do this for us::

    >>> my_datasets["RNA-seq"].set_mapping(CenterMapFactory())

Now, we need to create a place to hold our data. We'll use dictionary of lists.
The call to :func:`copy.deepcopy` on the empty list is necessary to prevent all
of these dictionary keys from pointing to the same list, which is a weird side
effect of the order in which things are evaluated inside comprehensions::

    >>> # we will count gene sub-regions in addition to entire genes
    >>> regions = ("exon","UTR5","CDS","UTR3")

    >>> # we will calculate both total counts and RPKM
    >>> # The `'%s'` notation is for string substitution. Each `%s` is substituted
    >>> # with the value of the variable in the tuple following it, in order from
    >>> # left-to-right. This is a convenient way to generate strings without using lots of
    >>> # string addition 
    >>> metrics = ("counts","rpkm")
    >>> my_data = { "%s_%s_%s" % (SAMPLE,REGION,METRIC) : copy.deepcopy([])\
    >>>                                                   for SAMPLE in datasets.keys()\
    >>>                                                   for REGION in regions\
    >>>                                                   for METRIC in metrics }

    >>> # add a list to our dictionary of lists to store transcript IDs
    >>> my_data["transcript_id"] = []

    >>> # add additional lists to store information about each region
    >>> for region in regions:
    >>>     my_data["%s_chain"  % region] = []  # SegmentChain representing region
    >>>     my_data["%s_length" % region] = []  # Length of that SegmentChain, in nucleotides


Now that we have an empty dictionary of lists to hold our data, we're ready to start
making measurements. We'll use nested for loops to count expression in the 5' UTR, 
CDS, 3'UTR and total region (exon) of each transcript::

    >>> for c, transcript in enumerate(GTF2_TranscriptAssembler(open("sgd_plus_utrs_chrI.gtf"))):
    >>> 
    >>>     if c % 100 == 0:
    >>>         print "Processed %s entries..." % c
    >>>     
    >>>     # First, save ID of transcript we are evaluating
    >>>     my_data["transcript_id"].append(transcript.get_name())

    >>>     # Next, get transcript sub-regions, save them in a dict
    >>>     # mapping region names to genomic regions (SegmentChains)
    >>>     my_dict = { "exon" : transcript,
    >>>                 "UTR5" : transcript.get_utr5(),
    >>>                 "CDS"  : transcript.get_cds(),
    >>>                 "UTR3" : transcript.get_utr3()
    >>>                }

    >>>     # Now, iterate over these sub-regions for each transcript
    >>>     for region,subchain in my_dict.items():
    >>>         # Save the length and a string representation of the coordinates
    >>>         # For each sub-region
    >>>         my_data["%s_length" % region].append(subchain.get_length())
    >>>         my_data["%s_chain"  % region].append(str(subchain))

    >>>         # Now, iterate over each sample, getting the counts under the 
    >>>         # mapping rules we set above. Here we iterate over the key-value
    >>>         # pairs in datasets, which map sample names to BAMGenomeArrays
    >>>         for sample_name, sample_data in datasets.items():
    >>>             # subchain.get_counts() fetches a list of counts at each position
    >>>             # here we just want the sum
    >>>             counts = sum(subchain.get_counts(sample_data))
    >>>             rpkm   = float(counts) / subchain.total_length() * 1000 * 1e6 / sample_data.sum()
    >>>             my_data["%s_%s_counts" % (sample_name,region)].append(counts)
    >>>             my_data["%s_%s_rpkm"   % (sample_name,region)].append(rpkm)

Finally, we can save the data to a file. It is easiest to do this by converting 
our dictionary of lists into a :class:`pandas.DataFrame`::

    >>> # convert to DataFrame, then save as tab-delimited text file
    >>> df = pd.DataFrame(my_data)
    >>> df.to_csv("%s_expression.txt" % sample,sep="\t")

That's it! These text files may be re-loaded for further analysis, or plotted.
For fun, let's plot the :term:`RPKM` measurements for translation
(:term:`ribosome profiling`) and transcription (:term:`RNA-seq`) against
each other::

    >>> my_figure = plt.figure()
    >>> plt.loglog() # log-scaling makes it easier

    >>> # make a copy of dataframe for plotting so we can set 0 values
    >>> # to a pseudo value (MIN_VAL), so they can be plotted in log space
    >>> MIN_VAL = 1e-5
    >>> plot_df = copy.deepcopy(df)
    >>> df["RNA-seq_exon_rpkm"][df["RNA-seq_exon_rpkm"] == 0] = MIN_VAL
    >>> df["ribosome_profiling_CDS_rpkm"][df["ribosome_profiling_CDS_rpkm"] == 0] = MIN_VAL


    >>> # now, make a scatter plot
    >>> plt.scatter(df["RNA-seq_exon_rpkm"],df["ribosome_profiling_CDS_rpkm"],
                    marker="o",alpha=0.2,facecolor="none",edgecolor="#007ADF")
    >>> plt.xlabel("Transcript levels (RPKM over all exons)")
    >>> plt.ylabel("Translation (RPKM of CDS)")

    >>> plt.show()



.. TODO : make image & insert here


Testing for differential expression
-----------------------------------

.. TODO : DESeq paper reference

There are many strategies for significance testing of differential gene expression
between multiple datasets. A very generalized and statistically rigorous approach
is taken by `DESeq`_, which takes as input the number of uncorrected :term:`counts`
in an arbitrary region of interest, in multiple datasets. 