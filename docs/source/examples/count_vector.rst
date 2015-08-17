Vectors of counts along transcripts
===================================
In this tutorial, we show:

  - :ref:`how to retrieve a vector <examples-count-vector-interactive>`
    of the number of 5' ends of :term:`read alignments` that align to each
    position in a transcript

  - :ref:`how to automate this process <examples-count-vector-script>`
    using the :mod:`~yeti.bin.get_count_vectors` script
 

In the examples below, we use the :doc:`/test_dataset`.


 .. _examples-count-vector-interactive:

Working with vectors of counts interactively
--------------------------------------------

 .. TODO : update all count vectors in this example

To count :term:`read alignments` along a transcript, we need two types of data:

  #. An :term:`annotation` of transcript models

  #. A :term:`high-throughput sequencing` dataset

First, we import everything we need::

    >>> # reader for BAM files
    >>> import pysam

    >>> # data structure for mapping read alignments to genomic positions
    >>> from yeti.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory

    >>> # reader for BED-format transcript annotations
    >>> from yeti.readers.bed import BED_Reader

Next, load the transcripts. By default, :class:`~yeti.genomics.readers.BED_Reader` 
and the other readers in :data:`yeti` behave as iterators. Here, we'll retrieve
the transcripts as a :class:`list`::

    >>> transcripts = list(BED_Reader(open("merlin_orfs.bed"),return_type=Transcript))

Then, load the :term:`ribosome profiling` data. The data are in a `BAM`_ file,
which we'll load into a :class:`~yeti.genomics.genome_array.BAMGenomeArray`.
We'll map :term:`read alignments` to the corresponding :term:`P-sites <P-site offset>`,
estimating the P-site to be 14 nucleotides from the 5' end::

    >>> alignments = BAMGenomeArray([pysam.Samfile("SRR609197_riboprofile.bam")])
    >>> alignments.set_mapping(FivePrimeMapFactory(offset=14))

Now, we're ready to count. The method
:meth:`Transcript.get_counts <yeti.genomics.roitools.Transcript.get_counts>`, returns
a vector (:class:`numpy.ndarray`) of counts corresponding to
each position in the transcript, from the 5' end of the transcript to the 3'
end (i.e. for reverse-strand features, counts are reversed relative to
genomic coordinates), accounting for splicing of exons::

    >>> count_vectors = []
    >>> for transcript in transcripts:
    >>>     count_vectors.append(transcript.get_counts(alignments))

    # we'll take transcript 53 as an example- it has lots of reads
    # check the lengths of the first transcript and its vector.
    # they should be identical
    >>> my_transcript = transcripts[53]
    >>> my_vector = count_vectors[53]
    >>> my_transcript.get_length(), len(my_vector)
    (1571, 1571)

    # get total counts over entire vector
    >>> my_vector.sum()
    7444.0

    >>> # slicing 
    >>> my_vector[200:250]
    array([   7.,   25.,   18.,   13.,    5.,    1.,   11.,    3.,    0.,
              1.,   25.,   11.,   29.,   27.,   18.,    3.,   16.,   20.,
             10.,    0.,    4.,   20.,   10.,    2.,    3.,   19.,    4.,
              9.,    1.,   15.,    5.,    3.,   11.,    8.,   13.,   15.,
              4.,  121.,    3.,    6.,   45.,    3.,    4.,   39.,   14.,
              3.,    9.,    7.,    8.,   24.])

Because the vector is a :class:`numpy.ndarray`, it can be manipulated using
any of the tools in `numpy`_, `SciPy`_, or `matplotlib`_::

    >>> import numpy
    
    # mean & variance in coverage
    >>> my_vector.mean(), my_vector.var()
    (4.7383831954169322, 49.177260021207104)

    # location of highest peak
    >>> my_vector.argmax()
    237

    # take cumulative sum
    >>> my_vector.cumsum()
    array([    0.,     0.,     0., ...,  7444.,  7444.,  7444.])
   
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
The analysis above is implemented by the command-line script |get_count_vectors|.
|get_count_vector| requires the same data types as above:

 #. An :term:`annotation` of genomic :term:`features <feature>`
    (e.g. transcripts for :term:`ribosome profiling`,
    promoters & enhancers for ChIP-seq, et c)
 
 #. Some :term:`high-throughput` sequencing data


The script may then be executed from the terminal:

 .. code-block:: shell

    $ get_count_vectors --annotation_files merlin_orfs.bed \
                        --annotation_format BED \
                        --count_files SRR609197_riboprofile.bam \
                        --fiveprime \
                        --offset 14 \
                        folder_of_vectors

Each output file will be saved in `folder_of_vectors` and named for the `ID`
attribute of the corresponding genomic :term:`feature`:

 .. code-block : shell                        

    $ ls folder_of_vectors
    ORFL100C.txt               ORFL169C.txt                 ORFL237C.txt                    ORFL308C_(UL139).txt         ORFL85C_(UL30).txt
    ORFL101C.iORF1_(UL36).txt  ORFL16C.iORF1.txt            ORFL238W.iORF1.txt              ORFL309C.txt                 ORFL86W.txt
    ORFL101C.txt               ORFL16C.txt                  ORFL238W.txt                    ORFL30W.txt                  ORFL87W.txt
    ORFL102C.iORF1.txt         ORFL170C.txt                 ORFL239C.txt                    ORFL310W.txt                 ORFL88C.iORF1.txt
    ORFL102C_(UL38).txt        ORFL171W.txt                 ORFL23W_(RL12).txt              ORFL311W.txt                 ORFL88C_(UL30A).txt
    ORFL103C_(vMIA).txt        ORFL172W.txt                 ORFL240C.txt                    ORFL312C.txt                 ORFL89C.txt
    ORFL104C_(UL37).txt        ORFL173W.txt                 ORFL241C_(UL103).txt            ORFL313C_(UL138).txt         ORFL8C.txt
    ORFL105C_(UL40).txt        ORFL174C.iORF2.txt           ORFL242W.txt                    ORFL314C.iORF1.txt           ORFL90C.txt
    (rest of output omitted)


The output can be loaded into numpy vectors using :func:`numpy.loadtxt`::

    >>> import numpy
    
    >>> my_reloaded_vector = numpy.loadtxt("folder_of_vectors/ORFL46W.iORF1_(UL13).txt")
    >>> my_reloaded_vector[200:250]
    array([   7.,   25.,   18.,   13.,    5.,    1.,   10.,    3.,    0.,
              1.,   24.,    9.,   27.,   27.,   18.,    3.,   16.,   20.,
             10.,    0.,    4.,   20.,   10.,    2.,    3.,   19.,    4.,
              9.,    1.,   15.,    5.,    3.,   11.,    8.,   13.,   14.,
              4.,  119.,    3.,    6.,   45.,    3.,    4.,   39.,   14.,
              3.,    9.,    7.,    8.,   24.])


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
