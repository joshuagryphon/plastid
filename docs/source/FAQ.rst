Frequently-asked questions
==========================

Questions are grouped into the following sections:

  - :ref:`faq-run`
  - :ref:`faq-analysis`
  - :ref:`faq-tests`


-------------------------------------------------------------------------------

 .. _faq-run:

Installation and runtime
------------------------

 .. _faq-install-numpy-first:

Install fails midway
....................

The most common cause of installation failure is that the setup script for
`plastid` requires `numpy`_, `pysam`_, and `Cython`_ to be installed before it
can be run. You can verify this by passing the ``--verbose`` flag to `pip`:

 .. code-block:: shell

    $ pip install --verbose plastid


You should see the following message:

 .. _faq-installation-message:

 .. code-block:: none

    *** IMPORTANT INSTALLATION INFORMATION ***

    plastid setup requires numpy>=1.9.0, pysam>=0.8.4, and cython>=0.22 to be preinstalled. Please
    install these via pip, and retry:

        $ pip install --upgrade numpy pysam cython
        $ pip install plastid

If this is the problem, simply install `numpy`_, `pysam`_ and `Cython`_ first,
then repeat the installation:

 .. code-block:: shell

    $ pip install numpy pysam cython # get these first
    $ pip install plastid


 .. _faq-install-versions:

`numpy`_, `pysam`_ and `Cython`_ are up to date, but install still fails
........................................................................

If install fails, first try the workaround in :ref:`faq-install-numpy-first`. If install
still fails with the message given in :ref:`faq-install-numpy-first`, then there may
be multiple versions of these libraries installed on your system. Try installing inside
a vanilla environment in a fresh `virtualenv`_ (note: *not* a `conda`_/`Anaconda`_ environment):

 .. code-block:: shell

    # install virtualenv if you don't have it
    $ pip install virtualenv

    # create & activate vanilla environment
    # when prompted, do NOT give the virtualenv access to system packages
    $ virtualenv /path/to/venv
    $ source path/to/venv/bin/activate

    # fresh install of plastid
    (venv) $ pip install numpy pysam cython
    (venv) $ pip install plastid

    # test
    (venv) $ python -c "from plastid import *"


If install succeeds, this suggests that there are in fact multiple versions of 
one or more of plastid's dependencies installed. In the mean time, ``plastid``
can be used inside the `virtualenv`_.


 .. _faq-install-conda:

Install fails inside a `conda`_/`Anaconda`_ environment
.......................................................

One user has reported difficulties installing inside `conda`_/`Anaconda`_
environments. Despite having up-to-date versions of `numpy`_, `pysam`_, and
`Cython`_ installed in the `conda`_ environment, during the build process
`pip` found an incompatible version of `Cython`_. In this case, we recommend
installing in a fresh `virtualenv`_. See the workaround in
:ref:`faq-install-versions` for instructions.


 .. _faq-distribution-error: 

I get an ``ImportError`` or ``DistributionError`` when using :data:`plastid`
............................................................................

If you get an error like the following::

    Traceback (most recent call last):
       File "/home/user/Rib_prof/venv/bin/crossmap", line 5, in <module>
         from pkg_resources import load_entry_point
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 2970, in <module>
         working_set = WorkingSet._build_master()
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 567, in _build_master
         ws.require(__requires__)
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 876, in require
         needed = self.resolve(parse_requirements(requirements))
       File "/home/user/Rib_prof/venv/lib/python2.7/site-packages/pkg_resources/__init__.py", line 761, in resolve
         raise DistributionNotFound(req)
     pkg_resources.DistributionNotFound: scipy>=0.12.0 


One or more dependencies (in this example, `SciPy`_ is not installed).
Please see :ref:`faq-install-numpy-first`, above.




-------------------------------------------------------------------------------

 .. _faq-analysis:
 
Analysis
--------


 .. _faq-analysis-fractional-counts:

Why do some scripts report fractional count numbers?
....................................................

Fractional counts for :term:`read alignments` arise when using a
alignment :term:`mapping rule` that maps reads fractionally over
multiple positions, for example to reflect uncertainty in the
exact position where the read should be counted. See the 
discussion of :doc:`concepts/mapping_rules`, where these are
discussed in depth.



 .. _faq-igv-vs-mapped-wiggle:

Why does `IGV`_ report way higher coverage at a given nucleotide than the file exported from |make_wiggle|?
...........................................................................................................
When `IGV`_ calculates coverage of a nucleotide, it counts the number of alignments covering that nucleotide.
So, a 30-nucleotide read would contribute 30 :term:`counts` to a dataset.

While it is possible to write any mapping rule in :mod:`plastid`, the :term:`mapping rules <mapping rule>`
included by default count each read only once (e.g. at their 5' end, 3' end, et c). Even when using
*center* or *entire* mapping, each position covered by a read alignment is only incremented by :math:`1.0/\ell`,
where :math:`\ell` is the length of the read. So, in this case, a 30-nucleotide read would only 
contribute 1 :term:`count <counts>` to a dataset. See :doc:`/concepts/mapping_rules/` for more information.


 .. _faq-cs-vs-counts-in-region:

What are the differences between :mod:`~plastid.bin.counts_in_region` and :mod:`~plastid.bin.cs`?
.................................................................................................
:mod:`~plastid.bin.counts_in_region` very simply counts read coverage (or any data) over
regions of interest, and reports those numbers in terms of :term:`counts` and :term:`RPKM`. It can 
optionally take a :term:`mask file`, if there are genomic positions in the regions
of interest which should be excluded from analysis. Otherwise, it makes no corrections.

:mod:`~plastid.bin.cs` is more complex, and is principally designed to make rough estimates
of gene expression at the gene, rather than transcript, level. In so doing, it makes several
heuristic corrections to regions before tabulating their :term:`counts` and :term:`RPKM`. Specifically:

 #. Genes that have transcripts that share exons are merged into single entities

 #. Gene areas are defined for each merged geen by including all positions occupied
    by all transcripts from that merged gene

 #. Regions occupied by two or more merged genes on the same strand are excluded from
    the calculation of expression values for both genes
 
 #. Optionally, a :term:`mask file` can be used to exclude any other positions from
    analysis.

 #. Expression values (in :term:`counts` and :term:`RPKM`) are tabulated for the entire
    gene area (reported as *exon_counts* and *exon_rpkm*) as well as for sub regions,
    if the gene is coding. Specifically, *cds_counts* and *cds_rpkm* are calculated
    from counts that cover positions in the gene area that are annotated as CDS in
    **all** transcripts in the merged gene. Ditto for 5' and 3' UTRs

Either one can be an appropriate starting place for a pipeline, depending upon your needs.
See the documentation and/or source code for |cs| and |counts_in_region| for further
discussion. 


 .. _faq-segmentchain-gff3:

Why does :meth:`plastid.genomics.roitools.SegmentChain.as_gff3` throw errors when exporting multi-segment chains?
..................................................................................................................
This is due to the incredible flexibility of the `GFF3`_ file format and ambiguities that this flexibility
necessarily induces. See :ref:`this advice <data-export-gff3>` on how to handle this.


 .. _faq-analysis-deseq:

How do I prepare data for differential gene expression analysis in `DESeq`_?
............................................................................

See :doc:`examples/gene_expression` in the :doc:`examples` section.


-------------------------------------------------------------------------------

 .. _faq-tests:

Tests
-----

The tests won't run
...................
In order to run the tests, you need to download the `test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_ and unpack it into ``plastid/test/``. We decided not to include the test data in the main package in order to keep the download small.




-------------------------------------------------------------------------------


 .. toctree::
    :maxdepth: 2
