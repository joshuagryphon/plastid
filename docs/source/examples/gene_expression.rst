Gene expression analysis
========================

This tutorial illustrates how to measure read density over regions. As 
an example, we look at gene expression (in :term:`raw read counts <counts>` and :term:`RPKM`)
using matched sampels of :term:`RNA-seq` and :term:`ribosome profiling` data.

However, the analysis below applies to other types of
:term:`high-throughput sequencing` data (e.g. ClIP-SEQ, :term:`DMS-Seq`, et c).

We will do this two ways:

 #. Using the :mod:`~yeti.bin.counts_in_region` script to
    :ref:`count expression automatically <gene-expression-scripts>`.

 #. :ref:`Manually <gene-expression-interactive>`, in an interactive
    Python session

The examples below use the :doc:`/test_dataset` we have assembled.


 .. note::

    For simplicity, we are performing these analyses without use 
  


Tabulating :term:`read counts <counts>` & :term:`RPKM`
------------------------------------------------------

 .. _gene-expression-scripts:

Via command-line scripts
........................

:data:`yeti` includes two scripts for measuring gene expression:

  * :mod:`~yeti.bin.cs`, which pre-processes a genome anntation and makes
    various heuristic corrections to gene boundaries (e.g. if genes overlap)

  * :mod:`~yeti.bin.counts_in_region`, which does not.

The differences between the scripts are further explained in
:ref:`faq-cs-vs-counts-in-region`. Here we will use :mod:`~yeti.bin.counts_in_region`.

Our first dataset is :term:`ribosome profiling`, and we will map the ribosomal
P-site at 14 nucleotides from the 5' end of each read (approximating :cite:`Stern-ginossar2012`).
To specify this, we use the arguments ``--fiveprime --offset 14``.

The data we want to count is in the file ``SRR609197_riboprofile.bam``, which we pass
via ``--count_files``. The genes we are interested in counting in this example
are on chromosome I, in the annotation file ``merlin_orfs.gtf``. Finally,
we will tell the script to save the output in ``riboprofile.txt``.

Putting this together, the script is run from the terminal as:

 .. code-block:: shell

    $ counts_in_region riboprofile.txt --count_files SRR609197_riboprofile.bam --annotation_files merlin_orfs.gtf --fiveprime --offset 14

:mod:`~yeti.bin.counts_in_region` will create a tab-delimited text file called
``riboprofile.txt`` containing the results. For detailed documentation of the output
and command-line arguments, see the module documentation for :mod:`~yeti.bin.counts_in_region`.

 .. note::

    When calculating densities of :term:`ribosome-protected footprints`
    over coding regions, the codons immediately surrounding the start
    and stop codons should be excluded from analysis due to the presence of
    :term:`start <start codon peak>` and :term:`stop codon peaks <stop codon peak>`,
    which inflate estimates of steady-state translation. We eliminated this
    step in this tutorial for the sake of simplicity.

    In practice, one could do this by supplying a :term:`mask annotation file <mask file>`
    annotating such regions. For a detailed explanation, see :doc:`/examples/using_masks`.


 .. _gene-expression-interactive:

Manually
........

Gene expression -- or, more broadly, read coverage over arbitrary regions
-- can be calculated easily in an interactive Python session.

In this example, we separately caclulate read density over:

  - entire transcripts
  - 5' UTRs
  - coding regions
  - 3' UTRs

First, we need to import a few things::

    >>> import copy
    >>> import pysam
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from yeti.readers.gff import GTF2_TranscriptAssembler
    >>> from yeti.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory, CenterMapFactory


Then, we'll open our data, storing each dataset in a |BAMGenomeArray|::

    >>> my_datasets = { "ribosome_profiling" : "SRR609197_riboprofile.bam",
    >>>                 "RNA-seq"            : "SRR592963_rnaseq.bam",
    >>>               }

    >>> my_datasets = { K : BAMGenomeArray([pysam.Samfile(V)]) for K,V in my_datasets.items() }

 
Next, we tell the |BAMGenomeArrays| which :term:`mapping rule` to use. We
will map the :term:`ribosome-protected footprints` to their P-sites, which
we estimate as 14 nucleotides from the 5' end of each read::

    >>> my_datasets["ribosome_profiling"].set_mapping(FivePrimeMapFactory(offset=14))

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
    >>> regions = ("exon","5UTR","CDS","3UTR")

    >>> # we will calculate both total counts and RPKM
    >>> metrics = ("counts","rpkm")

    >>> # create an empty list for each sample, region, and metric
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
CDS, 3'UTR and total region (exon) of each transcript:

 .. code-block:: python

    >>> for transcript in GTF2_TranscriptAssembler(open("merlin_orfs.gtf")):
    >>> 
    >>>     # First, save ID of transcript we are evaluating
    >>>     my_data["transcript_id"].append(transcript.get_name())

    >>>     # Next, get transcript sub-regions, save them in a dict
    >>>     # mapping region names to genomic regions (SegmentChains)
    >>>     my_dict = { "exon" : transcript,
    >>>                 "5UTR" : transcript.get_utr5(),
    >>>                 "CDS"  : transcript.get_cds(),
    >>>                 "3UTR" : transcript.get_utr3()
    >>>                }

    >>>     # Iterate over these sub-regions for each transcript
    >>>     for region,subchain in my_dict.items():
    >>>         # Save the length for each sub-region
    >>>         my_data["%s_length" % region].append(subchain.get_length())

    >>>         # Iterate over each sample, getting the counts over each region
    >>>         for sample_name, sample_data in datasets.items():
    >>>             # subchain.get_counts() fetches a list of counts at each position
    >>>             # here we just want the sum
    >>>             counts = sum(subchain.get_counts(sample_data))
    >>>             rpkm   = float(counts) / subchain.get_length() * 1000 * 1e6 / sample_data.sum()
    >>>             my_data["%s_%s_counts" % (sample_name,region)].append(counts)
    >>>             my_data["%s_%s_rpkm"   % (sample_name,region)].append(rpkm)


 .. note::

    As mentioned above, when calculating densities of :term:`ribosome-protected footprints`
    over coding regions, the codons immediately surrounding the start
    and stop codons should be excluded from analysis due to the presence of
    :term:`start <start codon peak>` and :term:`stop codon peaks <stop codon peak>`,
    which inflate estimates of steady-state translation.

    We skipped this here in order to keep this tutorial simple. For a detailed
    discussion of using masks in interactive Python sessions, see
    :ref:`masking-mask-file-interactive`.


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

    >>> # make a copy of dataframe for plotting
    >>> # this is because 0-values cannot be plotted in log-space,
    >>> # so we set them to a pseudo value called `MIN_VAL`
    >>>
    >>> MIN_VAL = 1e-5
    >>> plot_df = copy.deepcopy(df)
    >>> df["RNA-seq_exon_rpkm"][df["RNA-seq_exon_rpkm"] == 0] = MIN_VAL
    >>> df["ribosome_profiling_CDS_rpkm"][df["ribosome_profiling_CDS_rpkm"] == 0] = MIN_VAL


    >>> # now, make a scatter plot
    >>> plt.scatter(plot_df["RNA-seq_exon_rpkm"],
    >>>             plot_df["ribosome_profiling_CDS_rpkm"],
    >>>             marker="o",alpha=0.2,facecolor="none",edgecolor="#007ADF")
    >>> plt.xlabel("Transcript levels (RPKM of mRNA fragments over all exons)")
    >>> plt.ylabel("Translation (RPKM of footprints over CDS)")

    >>> plt.show()

This produces the following plot:

     .. figure:: 
        :figclass: captionfigure
        :alt: Scatter plot of translation versus transcription levels

        Translation versus transcription levels for each gene


Estimating translation efficiency
---------------------------------

:term:`Translation efficiency` is a measure of how much protein is
made from a single mRNA. :term:`Translation efficiency` thus reports
specifically on the *translational* control of gene expression.

:term:`translation efficiency` can be estimated
by normalizing an mRNA 's translating ribosome density (in :term:`RPKM`,
as measured by :term:`ribosome profiling`) by the mRNA's abundance (in
:term:`RPKM`, measured by :term:`RNA-Seq`) (:cite:`Ingolia2009`).

Making this estimate from the calculations above is trivial::

    >>> df["translation_efficiency"] = df["ribosome_profiling_CDS_rpkm"] / df["RNA-seq_exon_rpkm"]

Then, we can compare the effects of transcriptional and translational
control:

    >>> plot_df = copy.deepcopy(df)
    >>> df["RNA-seq_exon_rpkm"][df["RNA-seq_exon_rpkm"] == 0] = MIN_VAL
    >>> df["translation_efficiency"][df["translation_efficiency"] == 0] = MIN_VAL

    >>> # now, make a scatter plot
    >>> plt.scatter(plot_df["RNA-seq_exon_rpkm"],
    >>>             plot_df["translation_efficiency"],
    >>>             marker="o",alpha=0.2,facecolor="none",edgecolor="#007ADF")
    >>> plt.xlabel("Transcript levels (RPKM of mRNA fragments over all exons)")
    >>> plt.ylabel("Translation efficiency")

    >>> plt.show()

 .. TODO::

    Consider adding information about GTI-Seq or other TE estimates



Testing for differential expression
-----------------------------------

RNA-Seq
.......
There are many strategies for significance testing of differential gene expression
between multiple datasets, many of which are specifically developed for -- and
make statistical corrections that assume -- :term:`RNA-seq` data
(:cite:`some citation` ). In particuar, for :term:`RNA-seq` data, `cufflinks`_ and `kallisto`_
are fast and efficient, and don't require any preprocessing within :data:`yeti` at all.
For further information on using those packages, see their documentation.

Any :term:`high-throughput sequencing` experiment
.................................................
For other experimental data types -- e.g. :term:`ribosome profiling`, :term:`DMS-Seq`,
:term:`ChIP-Seq`, :term:`ClIP-Seq`, et c -- the assumptions made by many packages
developed for :term:`RNA-seq` analysis do not hold. 

In contrast, the `R`_ package `DESeq`_ offers a very generally applicable
statistical approach that is
appropriate to virtually any count-based sequencing data (:cite:`Anders2010`,
:cite:`Anders2013`).

As input, `DESeq` takes two objects:

 #. A table of *uncorrected, unnormalized* :term:`counts`, in which
    each row in the table corresponds to a genomic region, and each column an
    experimental sample. The value in each cell corresponds to the number of counts
    in that region, in that sample.

 #. An experimental design describing the relationships between samples
    (e.g. if any are technical or biological replicates)
     
From these, `DESeq`_ separately models intrinsic counting error as well
as inter-replicate error, and from these error models can infer significant
differences in count numbers between non-replicate samples.

The first table may be constructed by running |cs| or |counts_in_region|
on each biological sample, and extracting the relvant columns from their
output. One could do this in a spreadsheet program, in Python, or from
the terminal. For example:

 .. TODO: flesh out this example

 .. code-block:: shell

    $ counts_in_region ...
    $ counts_in_region ...
    $ counts_in_region ...
    $ counts_in_region ...

    $ # check which columns we want from each file
    $ head --line 20 ...
    
    $ echo "sample1\tsample2" >new_file.txt # note! replace '\t' with a tab when you type this!

    $ # create a new tab-delimited file from those columns
    $ paste <(cut -f ... ) <(cut -f ... ) >>new_file.txt


The second table is specified as an experimental *design*.

 .. code-block:: r


Differential translation efficiency
...................................

 .. TODO :

    Tests for differential translation efficiency can also be implemented within
    `DESeq`_..


Statistical models for differential measurement of :term:`translation efficiency`
are still a subject of discussion (TODO: citations). Here, we take an empirical
approach used in :cite:`Ingolia2009`.

 #. First, a :term:`false discovery rate` appropriate to the experiment 
    -- often five percent -- is set.

 #. For each sample, the :term:`translation efficiency` of each mRNA measured as
    the ratio of :term:`ribsome-protected footprint` density in a coding region
    to the mRNA fragment density across the corresponding mRNA.
 
 #. Within each set of biological replicates, log2 fold-changes are calculated
    for each transcript to yield an empirical distribution of changes derived
    from sequencing error for that replicate set. These distributions are 
    merged by summing the sets of their observations.

 #. Similarly, log2 fold-changes are calculated for each transcript between
    non-replicate samples. 

 #. The number of false positives (FP) at a given fold-change may be estimated
    as the number of observed fold changes greater to or equal than
    the given fold-change in the negative control distribution from step (3).

 #. Similarly, the number of total positives (FP+TP) at a given fold-change is the
    number of observed fold-changes greater to or equal than that fold-change
    in the distribution from step (4).

 #. The number of true positives (TP) at each fold-chnage is then estimated by subtracting
    the number of false positives at that fold-change (step 5) from the number
    of total positives (step 6).

 #. A significance threshold is set by solving for the fold change that corresponds
    to the :term:`false discovery rate (FDR) <false discovery rate>` set in step (1). 
    :term:`FDR` is calculated at each fold-change threshold :math:`t` as:

     .. math::

        FDR(t) = \frac{TP(t)}{TP(t)+FP(t)}

    Then, the fold-change :math:`t` where :term:`FDR` equals the predetermined
    :term:`false discovery rate` is taken to be the significance threshold.



-------------------------------------------------------------------------------

See also
--------

  - :doc:`/concepts/mapping_rules` and :mod:`yeti.genomics.genome_array` for
    information on mapping rules and processing read alignments

  - Documentation for |cs| and |counts_in_region| for further discussion 
    of their algorithms

  - `DESeq` website, :cite:`Anders2010`, and :cite:`Anders2013` for discussions
    on statistical models for differential gene expression, an examples
    on how to use `DESeq` for various experimental setups

  - :doc:`/examples/using_masks` for instructions on how to exclude parts of
    the genome or transcriptome from analysis.
