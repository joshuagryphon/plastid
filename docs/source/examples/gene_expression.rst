Gene expression analysis
========================

This tutorial illustrates how to measure read density over regions. As 
an example, we look at gene expression (in :term:`raw read counts <counts>` and :term:`RPKM`)
using matched samples of :term:`RNA-seq` and :term:`ribosome profiling` data.
However, the analysis below can apply to any type of
:term:`high-throughput sequencing` data (e.g. ClIP-seq, ChIP-seq, :term:`DMS-seq`, et c).

We will do this two ways:

 #. :ref:`Manually <gene-expression-interactive>`, in an interactive
    Python session

 #. Using the :mod:`~plastid.bin.counts_in_region` script to
    :ref:`count expression automatically <gene-expression-scripts>`.

The examples below the :doc:`/test_dataset`.


 .. note::

    For simplicity, we are:
    
     #. performing these analyses without using a :term:`mask file`
        to exclude repetitive genome seqence from analysis. In
        practice, we would consider doing so. See
        :doc:`/examples/using_masks` for a discussion of the
        issues involved

     #. not excluding start and stop codons from read density
        measurements in :term:`ribosome profiling` data. In practice,
        we would do so, because :term:`start <start codon peak>`
        and :term:`stop codon peaks <stop codon peak>`, artificially
        inflate estimates of steady-state translation.
        

Tabulating :term:`read counts <counts>` & :term:`RPKM`
------------------------------------------------------

 .. _gene-expression-scripts:

Via command-line scripts
........................

:data:`plastid` includes two scripts for measuring gene expression:

  * :mod:`~plastid.bin.cs`, which pre-processes a genome anntation and makes
    various heuristic corrections to gene boundaries (e.g. if genes overlap)

  * :mod:`~plastid.bin.counts_in_region`, which does not.

The differences between the scripts are further explained in
:ref:`faq-cs-vs-counts-in-region`. Here we will use :mod:`~plastid.bin.counts_in_region`.

Our first dataset is :term:`ribosome profiling`, and we will map the ribosomal
P-site at 14 nucleotides from the 5' end of each read (approximating :cite:`Stern-Ginossar2012`).
To specify this, we use the arguments ``--fiveprime --offset 14``.

The data we want to count is in the file ``SRR609197_riboprofile.bam``, which we pass
via ``--count_files``. The genes we are interested in counting in this example
are on chromosome I, in the annotation file ``merlin_orfs.gtf``. Finally,
we will tell the script to save the output in ``riboprofile.txt``.

Putting this together, the script is run from the terminal as:

 .. code-block:: shell

    $ counts_in_region riboprofile.txt --count_files SRR609197_riboprofile.bam \
                                       --annotation_files merlin_orfs.gtf \
                                       --fiveprime --offset 14

:mod:`~plastid.bin.counts_in_region` will create a tab-delimited text file called
``riboprofile.txt`` containing the results. The first few lines of the file
look like this::

    ## total_dataset_counts: 500477
    #region_name    region                  counts          counts_per_nucleotide   rpkm            length
    ORFL1W_(RL1)    merlin:1316-2398(+)     1.14000000e+02  1.05360444e-01          2.10520051e+02  1082
    ORFL2C          merlin:2401-2772(-)     1.00000000e+01  2.69541779e-02          5.38569762e+01  371
    ORFL3C          merlin:2834-3064(-)     1.50000000e+01  6.52173913e-02          1.30310466e+02  230
    ORFL4C          merlin:2929-3201(-)     1.40000000e+01  5.14705882e-02          1.02843064e+02  272
    ORFL5C          merlin:4074-4307(-)     2.30000000e+01  9.87124464e-02          1.97236729e+02  233
    ORFL6C          merlin:4078-4488(-)     6.10000000e+01  1.48780488e-01          2.97277373e+02  410
    ORFL7C          merlin:4335-4739(-)     6.20000000e+01  1.53465347e-01          3.06638160e+02  404
    [rest of output omitted]



For detailed documentation of the output and command-line arguments, see
the module documentation for :mod:`~plastid.bin.counts_in_region`.


 .. _gene-expression-interactive:

Manually
........

Gene expression -- or, more broadly, read density over from any
:term:`high-throughput sequencing` experiment over any genomic
region -- can be calculated easily in an interactive Python
session.

In this example, we separately caclulate read density over:

  - entire transcripts
  - 5' UTRs
  - coding regions
  - 3' UTRs

First, we need to import a few things::

    >>> import copy
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> from plastid import Transcript BED_Reader, \
                            BAMGenomeArray, FivePrimeMapFactory, CenterMapFactory


First, open the :term:`read alignments`, storing each dataset in a |BAMGenomeArray|::

    >>> my_datasets = { "ribosome_profiling" : "SRR609197_riboprofile.bam",
    >>>                 "RNA-seq"            : "SRR592963_rnaseq.bam",
    >>>               }

    >>> my_datasets = { K : BAMGenomeArray([V]) for K,V in my_datasets.items() }

 
Next, we tell the |BAMGenomeArrays| which :term:`mapping rule` to use. We
will map the :term:`ribosome-protected footprints` to their P-sites, which
we estimate as 14 nucleotides from the 5' end of each read::

    >>> my_datasets["ribosome_profiling"].set_mapping(FivePrimeMapFactory(offset=14))

We will map the RNA-seq data along the entire length of each read alignment.
Each position in each alignment will be attributed :math:`1.0 / \ell`, where 
:math:`\ell` is the length of the read alignment.
:func:`~plastid.genomics.genome_array.CenterMapFactory` can do this for us::

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
    >>>                                                   for SAMPLE in my_datasets.keys()\
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
CDS, 3'UTR and total region (exon) of each transcript (**note:** this will run for a 
while; you might want to get some coffee):

 .. code-block:: python

    >>> for transcript in BED_Reader(open("merlin_orfs.bed"),return_type=Transcript):
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
    >>>         my_data["%s_length" % region].append(subchain.length)
    >>>         my_data["%s_chain"  % region].append(str(subchain))

    >>>         # Iterate over each sample, getting the counts over each region
    >>>         for sample_name, sample_data in my_datasets.items():
    >>>             # subchain.get_counts() fetches a list of counts at each position
    >>>             # here we just want the sum
    >>>             counts = sum(subchain.get_counts(sample_data))
    >>>             rpkm   = float(counts) / subchain.length * 1000 * 1e6 / sample_data.sum()
    >>>             my_data["%s_%s_counts" % (sample_name,region)].append(counts)
    >>>             my_data["%s_%s_rpkm"   % (sample_name,region)].append(rpkm)

Finally, we can save the calculated values to a file. It is easiest to do this
by converting the dictionary of lists into a :class:`pandas.DataFrame`:: 

    >>> # convert to DataFrame, then save as tab-delimited text file
    >>> df = pd.DataFrame(my_data)
    >>> df.to_csv("gene_expression_demo.txt",sep="\t",index=False,header=True)

The text files may be re-loaded for further analysis, or plotted. For example,
to plot the :term:`RPKM` measurements for translation (:term:`ribosome profiling`)
and transcription (:term:`RNA-seq`) against each other::

    >>> my_figure = plt.figure()
    >>> plt.loglog() # log-scaling makes it easier

    >>> # make a copy of dataframe for plotting
    >>> # this is because 0-values cannot be plotted in log-space,
    >>> # so we set them to a pseudo value called `MIN_VAL`
    >>>
    >>> MIN_VAL = 1
    >>> plot_df = copy.deepcopy(df)
    >>> df["RNA-seq_exon_rpkm"][df["RNA-seq_exon_rpkm"] == 0] = MIN_VAL
    >>> df["ribosome_profiling_CDS_rpkm"][df["ribosome_profiling_CDS_rpkm"] == 0] = MIN_VAL

    >>> # now, make a scatter plot
    >>> plt.scatter(plot_df["RNA-seq_exon_rpkm"],
    >>>             plot_df["ribosome_profiling_CDS_rpkm"],
    >>>             marker="o",alpha=0.5,facecolor="none",edgecolor="#007ADF")
    >>> plt.xlabel("Transcript levels (RPKM of mRNA fragments over all exons)")
    >>> plt.ylabel("Translation (RPKM of footprints over CDS)")

    >>> plt.show()


This produces the following plot:

     .. figure:: /_static/images/demo_gene_expr_tl_vs_tx.png
        :figclass: captionfigure
        :alt: Scatter plot of translation versus transcription levels

        Translation versus transcription levels for each gene


Estimating translation efficiency
---------------------------------

:term:`Translation efficiency` is a measurement of how much protein is
made from a single mRNA. :term:`Translation efficiency` thus reports
specifically on the *translational* control of gene expression.

:term:`Translation efficiency` can be estimated
by normalizing an mRNA 's translating ribosome density (in :term:`RPKM`,
as measured by :term:`ribosome profiling`) by the mRNA's abundance (in
:term:`RPKM`, measured by :term:`RNA-Seq`) (:cite:`Ingolia2009`).

Making this estimate from the calculations above is simple::

    >>> df["translation_efficiency"] = df["ribosome_profiling_CDS_rpkm"] / df["RNA-seq_exon_rpkm"]

Then, we can compare the effects of transcriptional and translational
control:

 .. code-block:: python

    >>> plt.loglog()
    >>> plot_df = copy.deepcopy(df)
    >>> plot_df["RNA-seq_exon_rpkm"][df["RNA-seq_exon_rpkm"] == 0] = MIN_VAL
    >>> plot_df["translation_efficiency"][df["translation_efficiency"] == 0] = MIN_VAL

    >>> # now, make a scatter plot
    >>> plt.scatter(plot_df["RNA-seq_exon_rpkm"],
    >>>             plot_df["translation_efficiency"],
    >>>             marker="o",alpha=0.2,facecolor="none",edgecolor="#007ADF")
    >>> plt.xlabel("Transcript levels (RPKM of mRNA fragments over all exons)")
    >>> plt.ylabel("Translation efficiency")
    >>> plt.xlim(1,plt.get_xlim()[1])
    >>> plt.ylim(plt.ylim()[0]/10.0,100)

    >>> plt.show()





 .. figure:: /_static/images/demo_gene_expr_teff_vs_tx.png

    :class: captionfigure
    :caption: Translation efficiency vs transcription levels
    :alt: Translation efficiency vs transcription levels



Testing for differential expression
-----------------------------------

 .. note::

    We need to add some more samples to the test dataset for this section.
    In the mean time, just read along


RNA-seq, specifically
.....................
There are many strategies for significance testing of differential gene expression
between multiple datasets, many of which are specifically developed for -- and
make statistical corrections that assume -- :term:`RNA-seq` data.

For :term:`RNA-seq` data, `cufflinks`_ and `kallisto`_ in particular are popular,
and operate directly on alignments in `BAM`_ format. These packages don't require
:data:`plastid` at all. For further information on them packages, see their documentation.


Any :term:`high-throughput sequencing` experiment, including RNA-seq
....................................................................
For other experimental data types -- e.g. :term:`ribosome profiling`, :term:`DMS-seq`,
:term:`ChIP-Seq`, :term:`ClIP-Seq`, et c -- the assumptions made by many packages
specifically developed for :term:`RNA-seq` analysis do not hold. 

In contrast, the `R`_ packages `DESeq`_ and `DESeq2`_ (:cite:`Anders2010,Anders2013,Love2014`)
offer a generally applicable statistical approach that is appropriate to virtually
any count-based sequencing data.

 .. note::
 
    The discussion below is heavily simplified and largely draws upon guidance in
    `Analysing RNA-Seq data with the "DESeq2" package <http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>`_,
    hosted on the `DESeq2`_ website.
    
    Users are encouraged to read the `DESeq`_/`DESeq2`_ documentation for a fuller
    discussion with additional examples.

As input, `DESeq`_ and `DESeq2`_ take two tables and an equation:

 #. A :ref:`table <examples-deseq-count-table>` of *uncorrected, unnormalized*
    :term:`counts`, in which:

      - each table row corresponds to a genomic region
      - each column corresponds to an experimental sample
      - the value in a each cell corresponds ot the number of counts
        in the corresponding genomic region and sample

 #. A :ref:`sample design table <examples-deseq-design-table>`
    describing the properties of each sample
    (e.g. if any are technical or biological replicates, or any treatments
    or conditions that differ between samples)

 #. A :ref:`design equation <examples-deseq-equation>`, describing how
    the samples or treatments relate to one another

    
From these, `DESeq`_ and `DESeq2`_ separately model intrinsic counting
error (Poisson noise) as well as additional inter-replicate error
resulting biological or experimental variability. From these error models,
`DESeq`_ and `DESeq2`_ can detect significant differences in count numbers
between non-replicate samples, accounting for different sequencing depth
between samples.


 .. _examples-deseq-count-table

The first table may be constructed by running |cs| or |counts_in_region|
on each biological sample to obtain counts. Here, we'll use RNA-seq data:

 .. code-block:: shell

    $ counts_in_region rnaseq_5hr_rep1 --count_files SRR592963_rnaseq_5hr_rep1.bam        --fiveprime --annotation_files merlin_orfs.bed --annotation_format BED
    $ counts_in_region rnaseq_5hr_rep2 --count_files _rnaseq_5hr_rep2.bam                 --fiveprime --annotation_files merlin_orfs.bed --annotation_format BED


    $ counts_in_region rnaseq_24hr_rep1 --count_files SRR592964_rnaseq_24hr_rep1.bam      --fiveprime --annotation_files merlin_orfs.bed --annotation_format BED
    $ counts_in_region rnaseq_24hr_rep2 --count_files SRR592967_rnaseq_24hr_rep2.bam --fiveprime --annotation_files merlin_orfs.bed --annotation_format BED



 .. TODO: include output
From the output, the relevant columns can be extracted and moved to
a single table:

 .. code-block:: python

    >>> import pandas as pd
    >>> sample_names = ["rnaseq_5hr_rep1","rnaseq_5hr_rep2","rnaseq_24hr_rep1","rnaseq_24hr_rep2"]

    >>> # load samples as DataFrames
    >>> samples = { K : pd.read_table("%s.txt" % K,sep="\t",header=0,comment="#",index_col=None) for K in sample_names }

    >>> # combine count columns to single DataFrame
    >>> combined_df = samples["rnaseq_5hr_rep1"]["region_name","region"]
    >>> for k,v in samples.items():
    >>>     combined_df["%s_counts" % k] = v["counts"]

    >>> combined_df.head()

    >>> # save
    >>> combined_df.savecsv("combined_counts.txt",sep="\t",header=True,index=False)






 .. _examples-deseq-design-table:

The second table contains the *experimental design*. This can be created
in any text editor and saved as a tab-delimited text file. In this example,
the we have two conditions, "infected" and "uninfected", and two replicates
of each condition:

 .. code-block :: shell

    sample_name          condition
    rnaseq_5hr_rep1      5_hours
    rnaseq_5hr_rep2      5_hours
    rnaseq_24hr_rep1     24_hours
    rnaseq_24hr_rep2     24_hours





 .. _examples-deseq-equation:

Because the only difference between samples is the `condition` column,
the design equation is this case is very simple::

    design = ~ condition


With the count table, design table, and equation ready, everything can
be loaded into `R`_:

 .. TODO: put output below
 .. code-block:: r

    > # load RNA seq data into a data.frame
    > # first line of file are colum headers
    > # "region" column specifies a list of row names
    > count_table <- read.delim("combined_counts.txt",
    >                           sep="\t",
    >                           header=True,
    >                           row.names="region")

    > sample_table <- read.delim("rnaseq_sample_table.txt",
    >                            sep="\t",
    >                            header=True,
    >                            row.names="sample_name")

    > # import DESeq2 & run with default settings
    > library("DESeq2")

    > # note, design string below tells DESeq2 that the 'condition' column
    > # distinguishes replicates from non-replicates 
    > dds <- DESeqDataSetFromMatrix(countData = count_table,
    >                               colData = sample_table,
    >                               design = ~ condition) # <--- design equation

    > results <- results(dds)
    > summary(res)

    > # sort results by adjusted p-value
    > resOrdered <- res[order(res$padj),]

    > # export sorted data to text file
    > write.delim(as.data.frame(resOrdered),
    >             sep="\t",
    >             file="infected_uninfected_rnaseq_p_values.txt")


Differential translation efficiency
...................................

Tests for differential translation efficiency can also be implemented within
`DESeq`_/`DESeq2`_. The discussion below follows a reply from `DESeq2`_ author
Mike Love (source `here <https://support.bioconductor.org/p/56736/>`_).

We use an sample table similar to that above, but include a `sample_type`
column to distinguish :term:`ribosome profiling` from :term:`RNA-seq` libraries::

    sample_name          condition    sample_type
    rnaseq_5hr_rep1      5_hours      rnaseq
    rnaseq_5hr_rep2      5_hours      rnaseq
    rnaseq_24hr_rep1     24_hours     rnaseq
    rnaseq_24hr_rep2     24_hours     rnaseq
    ribo_5hr_rep1        5_hours      riboprof
    ribo_24hr_rep1       24_hours     riboprof


To the design equation, we need to add  an *interaction term* to alert
`DESeq`_/`DESeq2`_ that we expect the relationship between the sample
types (i.e. translation efficiency, the ratio of
:term:`ribosome-protected footprints <footprint>` to RNA-seq fragments)
to differ between conditions::

    design = ~ sample_type + condition + sample_type:condition

In `R`_:

 .. TODO: put output below
 .. code-block:: r

    > # load RNA seq data into a data.frame
    > # first line of file are colum headers
    > # "region" column specifies a list of row names
    > combined_data <- read.delim("combined_counts.txt",
    >                             sep="\t",
    >                             header=True,
    >                            row.names="region")

    > teff_sample_table <- read.delim("teff_sample_table.txt",
    >                                sep="\t",
    >                                header=True,
    >                                row.names="sample_name")

    > library("DESeq2")

    > # note the interaction term in the design below:
    > dds <- DESeqDataSetFromMatrix(countData = combined_data,
    >                               colData = teff_sample_table,
    >                               design = ~ sample_type + condition + sample_type:condition)

    > results <- results(dds)
    > summary(res)

    > # now, do wald test on interaction term
    TODO: complete this line

    > # sort by adjusted p-value
    > resOrdered <- res[order(res$padj),]

    > # export
    > write.delim(as.data.frame(resOrdered),
    >             sep="\t",
    >             file="infected_uninfected_rnaseq_p_values.txt")


 .. old discussion- the empirical test used by Nick Ingolia 

    Statistical models for differential measurement of :term:`translation efficiency`
    are still a subject of discussion (TODO: citations). Here, we take an empirical
    approach used in :cite:`Ingolia2009`.

     #. First, a :term:`false discovery rate` (:cite:`Benjamini1995`) appropriate
        to the experiment -- often five percent -- is set.

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

  - :doc:`/concepts/mapping_rules` and :mod:`plastid.genomics.genome_array` for
    information on mapping rules and processing read alignments

  - Documentation for |cs| and |counts_in_region| for further discussion 
    of their algorithms

  - Websites for `DESeq` and `DESeq2`_, as well as :cite:`Anders2010`,
    :cite:`Anders2013` and :cite:`Love2014` for discussions of statistical models
    for differential gene expression, an examples
    on how to use `DESeq`_/`DESeq2`_ for various experimental setups

  - :doc:`/examples/using_masks` for instructions on how to exclude parts of
    the genome or transcriptome from analysis.
