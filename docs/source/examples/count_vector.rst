Vectors of counts along transcripts
===================================
In this tutorial, we show:

  - :ref:`how to retrieve a vector <examples-count-vector-interactive>`
    of the number of :term:`counts` from a
    :term:`high-throughput sequencing` experiment that align to each
    position in a transcript, in an interactive Python session

  - :ref:`how to automate this process <examples-count-vector-script>`
    using the :mod:`~yeti.bin.get_count_vectors` script
 

In the examples below, we use the :doc:`/test_dataset`.


 .. _examples-count-vector-interactive:

Working with vectors of counts interactively
--------------------------------------------

This tutorial requires two types of data:

  #. An :term:`annotation` of transcript models.

  #. A :term:`high-throughput sequencing` dataset.

First, we import everything we need::

    >>> import pysam 
    >>> from yeti.genomics.genome_array import BAMGenomeArray, ThreePrimeMapFactory
    >>> from yeti.readers.bed import BED_Reader

Next, load the transcripts. By default, :class:`~yeti.genomics.readers.BED_Reader` 
and the other readers in :data:`yeti` behave as iterators. Here, we'll retrieve
the transcripts as a :class:`list`::

    >>> reader = list(BED_Reader(open("sgd_plus_utrs_chrI.bed")))

Then, load the :term:`ribosome profiling` data. The data are in a `BAM`_ file,
which we'll load into a :class:`~yeti.genomics.genome_array.BAMGenomeArray`.
For our :term:`mapping rule`, we'll map reads 15 nucleotides from the 3' end,
corresponding to the P-site of the ribosome (cite:`Ingolia2009`)::

    >>> alignments = BAMGenomeArray([pysam.Samfile("SRR1562907_chrI.bam")])
    >>> alignments.set_mapping(ThreePrimeMapFactory(offset=15))

Then, to fetch a vector of counts covering each transcript, we'll use
the :meth:`Transcript.get_counts <yeti.genomics.roitools.Transcript.get_counts>`,
which returns a fully-spliced vector (:class:`numpy.ndarray`) of counts corresponding to
each position in the transcript, from the 5' end of the transcript to the 3'
end (i.e. for reverse-strand features, counts are reversed relative to
genomic coordinates)::

    >>> count_vectors = []
    >>> for transcript in transcripts:
    >>>     count_vectors.append(transcript.get_counts(alignments))

    # we'll take transcript 38 as an example- it has lots of reads
    # check the lengths of the first transcript and its vector.
    # they should be identical
    >>> my_transcript = transcripts[38]
    >>> my_vector = count_vectors[38]
    >>> my_transcript.get_length(), len(my_vector)
    (1698, 1698)

    # get total counts over entire vector
    >>> my_vector.sum()
    367301.0

    >>> # slicing - look at first 100 positions of vector
    >>> my_vector[:100]
    array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             2.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   6.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   1.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   1.00000000e+00,   2.00000000e+00,
             3.00000000e+00,   1.00000000e+00,   0.00000000e+00,
             2.00000000e+00,   0.00000000e+00,   1.40000000e+01,
             2.30000000e+01,   1.40000000e+02,   1.78500000e+03,
             3.34000000e+02,   4.40000000e+01,   7.53000000e+02,
             4.10000000e+02,   1.20000000e+02,   2.49000000e+02,
             1.11000000e+02,   2.99000000e+02,   1.45300000e+03,
             1.76000000e+02,   4.10000000e+01,   2.46000000e+02,
             5.99000000e+02,   4.80000000e+01,   5.59000000e+02,
             7.60000000e+01,   5.90000000e+01,   8.73000000e+02,
             3.70000000e+01,   4.50000000e+01,   6.48000000e+02,
             2.50000000e+01,   5.60000000e+01,   6.23000000e+02,
             4.60000000e+01,   3.60000000e+01,   5.63000000e+02,
             6.90000000e+01,   1.10000000e+01,   8.45000000e+02,
             1.38000000e+02,   6.80000000e+01,   6.25000000e+02,
             2.20000000e+01,   1.50000000e+01,   3.81000000e+02,
             2.35000000e+02,   6.50000000e+01,   4.07000000e+02,
             2.26000000e+02,   9.50000000e+01,   6.75000000e+02,
             1.59000000e+02,   5.40000000e+01,   8.46000000e+02,
             8.50000000e+01,   9.50000000e+01,   4.33000000e+02,
             1.06000000e+02,   1.23000000e+02,   5.97000000e+02,
             4.10000000e+01,   1.91000000e+02,   2.83000000e+02,
             8.60000000e+01,   7.20000000e+01,   4.20000000e+02,
             7.40000000e+01,   3.10000000e+01,   6.86000000e+02,
             4.30000000e+01,   7.70000000e+01,   4.60000000e+02,
             4.60000000e+01,   2.80000000e+01,   1.47000000e+02,
             1.47000000e+02,   3.40000000e+01,   7.77000000e+02,
             1.02000000e+02])

Because the vector is a :class:`numpy.ndarray`, it can be manipulated using
any of the tools in `numpy`_, `SciPy`_, or `matplotlib`_::

    >>> import numpy
    
    # mean & variance in coverage
    >>> my_vector.mean(), my_vector.var()
    >>> (216.31389870435808, 86983.56872319823)

    # location of highest peak
    >>> my_vector.argmax()
    218

    # take cumulative sum
    >>> my_vector.cumsum()
    array([      0.,       0.,       0., ...,  367301.,  367301.,  367301.])
   
    # 30-codon sliding window average
    >>> window = numpy.ones(90).astype(float)/90.0
    >>> sliding_window_avg = numpy.convolve(my_vector,window,mode="valid")


    # plot
    >>> import matplotlib.pyplot as plt

    >>> plt.plot(my_vector,label="%s counts" % my_transcript.get_name())
    >>> plt.plot(sliding_window_avg,label="30 codon average")
    >>> plt.xlabel("Position in transcript (5' to 3')")
    >>> plt.ylabel("Ribosome counts")

    >>> # add outlines at start & stop codons
    >>> plt.axvline(my_transcript.cds_start,color="#999999",dashes=[3,2],zorder=-1)
    >>> plt.axvline(my_transcript.cds_end,color="#999999",dashes=[3,2],zorder=-1)

    >>> plt.legend()
    >>> plt.show()

This makes the following figure:

 .. figure:: /_static/images/count_vectors_transcript_plot.png
    :figclass: captionfigure
    :alt: Sample plot of ribosome density

    Ribosome density at each position in a sample transcript. Dashed vertical lines:
    start and stop codons.


 .. _examples-count-vector-script:

Using the |get_count_vectors| script
------------------------------------
The analysis above is performed by the command-line script
|get_count_vectors|.

To run, this script requires the same
data types as above:

 #. An :term:`annotation` of genomic :term:`features <feature>`
    (e.g. transcripts for :term:`ribosome profiling`,
    promoters & enhancers for ChIP-seq, et c)
 
 #. Some :term:`high-throughput` sequencing data


The script may then be executed from the terminal:

 .. code-block:: shell

    $ get_count_vectors --annotation_files sgd_plus_utrs_chrI.bed \
                        --annotation_format BED \
                        --count_files SRR1562907_chrI.bam \
                        --threeprime --offset 15 \
                        folder_of_vectors

Each output file will be saved in `folder_of_vectors` and named for the `ID`
attribute of the corresponding genomic :term:`feature`:

 .. code-block : shell                        
    $ ls folder_of_vectors
    HRA1.txt            YAL019W_mRNA.txt    YAL039C_mRNA.txt    YAL063C-A_mRNA.txt  YAR027W_mRNA.txt
    snR18.txt           YAL020C_mRNA.txt    YAL040C_mRNA.txt    YAL063C_mRNA.txt    YAR028W_mRNA.txt
    tA(UGC)A.txt        YAL021C_mRNA.txt    YAL041W_mRNA.txt    YAL064C-A_mRNA.txt  YAR029W_mRNA.txt
    tL(CAA)A.txt        YAL022C_mRNA.txt    YAL042C-A_mRNA.txt  YAL064W-B_mRNA.txt  YAR030C_mRNA.txt
    tP(UGG)A.txt        YAL023C_mRNA.txt    YAL042W_mRNA.txt    YAL064W_mRNA.txt    YAR031W_mRNA.txt
    tS(AGA)A.txt        YAL024C_mRNA.txt    YAL043C_mRNA.txt    YAL065C_mRNA.txt    YAR033W_mRNA.txt
    YAL001C_mRNA.txt    YAL025C_mRNA.txt    YAL044C_mRNA.txt    YAL066W_mRNA.txt    YAR035C-A_mRNA.txt
    (rest of output omitted)    


The output can be loaded into numpy vectors using :func:`numpy.loadtxt`::

    >>> import numpy
    
    >>> my_reloaded_vector = numpy.loadtxt("folder_of_vectors/YAL038W_mRNA.txt")
    >>> my_reloaded_vector[:30]
    array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             2.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   6.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   1.00000000e+00,   0.00000000e+00,
             1.00000000e+00,   1.00000000e+00,   2.00000000e+00,
             3.00000000e+00,   1.00000000e+00,   0.00000000e+00,
             2.00000000e+00,   0.00000000e+00,   1.40000000e+01,
             2.30000000e+01,   1.40000000e+02,   1.78500000e+03])


|get_count_vectors| can optionally take a :term:`mask file` to exclude
problematic regions from analysis. In this case, vectors are returned
as :class:`numpy.ma.MaskedArray` objects, and positions annotated
in the :term:`mask file` are given the value :obj:`numpy.NaN` instead
of their numerical values. See :doc:`/examples/using_masks` for a 
discussion of :term:`mask files <mask file>` and how to make them
using |crossmap|.

-------------------------------------------------------------------------------

See also
--------
  - :doc:`/concepts/mapping_rules` for further discussion of
    :term:`mapping rules <mapping rule>`

  - :class:`~yeti.genomics.genome_array.GenomeArray` and
    :class:`~yeti.genomics.genome_array.BAMGenomeArray` for
    descriptions of Genome Arrays

  - :class:`~yeti.genomics.roitools.SegmentChain` and
    :class:`~yeti.genomics.roitools.Transcript` for full documentation
    of what these objects can do

  - :mod:`yeti.readers` subpackage, for readers
    of other :term:`annotation` file formats
