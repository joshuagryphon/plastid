Count vectors along transcripts
===============================
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

    import pysam 
    from yeti.genomics.genome_array import BAMGenomeArray, ThreePrimeMapFactory
    from yeti.genomics.readers.bed import BED_Reader

Next, load the transcripts. By default, :class:`~yeti.genomics.readers.BED_Reader` 
and the other readers in :data:`yeti` behave as iterators. Here, we'll retrieve
the transcripts as a :class:`list`::

    reader = list(BED_Reader(open("sgd_plus_utrs_chrI.bed")))

Then, load the :term:`ribosome profiling` data. The data are in a `BAM`_ file,
which we'll load into a :class:`~yeti.genomics.genome_array.BAMGenomeArray`.
For our :term:`mapping rule`, we'll map reads 15 nucleotides from the 3' end,
corresponding to the P-site of the ribosome (cite:`Ingolia2009`)::

    alignments = BAMGenomeArray([pysam.Samfile("SRR1562907_chrI.bam")])
    alignments.set_mapping(ThreePrimeMapFactory(offset=15))

Then, to fetch a vector of counts covering each transcript, we'll use
the :meth:`Transcript.get_counts <yeti.genomics.roitools.Transcript.get_counts>`,
which returns a fully-spliced vector (:class:`numpy.ndarray`) of counts corresponding to
each position in the transcript, from the 5' end of the transcript to the 3'
end (i.e. for reverse-strand features, counts are reversed relative to
genomic coordinates)::

    count_vectors = []
    for transcript in transcripts:
        count_vectors.append(transcript.get_counts(alignments))

    # check the lengths of the first transcript and its vector.
    # they should be identical
    transcripts[0].get_length(), len(count_vectors[0])

    # get total counts over entire vector
    # count_vectors[0].sum()

    # look at vector
    count_vectors[0]


Because the vector is a :class:`numpy.ndarray`, it can be manipulated using
any of the tools in :data:`numpy`, :data:`scipy`, or :data:`matplotlib`::

    my_transcript = transcripts[0]
    my_vector = count_vectors[0]

    # mean & variance
    my_vector.mean(), my_vector.var()

    # location of highest peak
    my_vector.argmax()


    # take cumulative sum
    my_vector.cumsum()

    # slice
    my_vector[200:300]


    # plot
    import matplotlib.pyplot as plt

    plt.plot(my_vector,label=my_transcript.get_name)
    plt.xlabel("Position in transcript (5' to 3')")
    plt.ylabel("Ribosome density")

    # add outlines at start & stop codons
    plt.axvline(my_transcript.cds_start,color="#999999",dashes=[3,2],zorder=-1)
    plt.axvline(my_transcript.cds_stop,color="#999999",dashes=[3,2],zorder=-1)

    plt.legend()

This makes the following figure:

 .. TODO: test examples
 .. TODO: make figure

 .. figure:: /_static/images/count_vectors_transcript_plot.png
    :figclass: captionfigure
    :alt: Sample plot of ribosome density

    Ribosome density at each position in a sample transcript. Dashed vertical lines:
    start and stop codons.



 .. _examples-count-vector-script:

Using the |get_count_vectors| script
------------------------------------------------------------
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

    $ get_count_vectors --annotation_files sgd_plus_utrs_chrI.gtf \
                        --annotation_format BED \
                        --count_files SRR1562907_chrI.bam \
                        --threeprime --offset 15 \
                        folder_of_vectors


|get_count_vectors| creates a folder, and saves each feature's
vector as a separate file in that folder. The output can be loaded
into numpy vectors using :func:`numpy.loadtxt`::

    import numpy
    
    my_vector = numpy.loadtxt("folder_of_vectors/some_vector.txt")
    my_vector


|get_count_vectors| can optionally take a :term:`mask file` to exclude
problematic regions from analysis. In this case, vectors are returned
as :class:`numpy.ma.masked_array` objects, and positions annotated
in the :term:`mask file` are given the value :obj:`numpy.nan` instead
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

  - :mod:`~yeti.genomics.readers` subpackage, for readers
    of other :term:`annotation` file formats
