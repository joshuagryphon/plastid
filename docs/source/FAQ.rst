Frequently-asked questions
==========================


-------------------------------------------------------------------------------

.. _faq-run:

Installation and runtime
------------------------

.. _faq-install-numpy-first:


The tests won't run
...................
In order to run the tests, you need to download the `test dataset <https://www.dropbox.com/s/h17go7tnas4hpby/plastid_test_data.tar.bz2?dl=0>`_ and unpack it into ``plastid/test/``. We didn't include this 
in the main package in order to keep the download small.


Install fails midway
....................

This is a due to a `known bug <https://github.com/numpy/numpy/issues/2434>`_ 
with `NumPy`_, `SciPy`_, and `matplotlib`_ installation under setuptools.
As a workaround. Just install `NumPy`_ first. Starting in the folder
where you cloned :data:`plastid`:

.. code-block:: shell

    $ pip install numpy cython # get these first
    $ pip install -r requirements.txt  # get everything else

You can your installations by running unit tests:

.. code-block:: shell

    (venv)$ python -c "import numpy; numpy.test()"
    (venv)$ python -c "import scipy; scipy.test()"


And make sure other packages are just importable:

.. code-block:: shell

    (venv)$ python -c "import matplotlib"
    (venv)$ python -c "import pysam"


Then repeat the installation:

.. code-block:: shell

    $ python setup.py build_ext --inplace # build plastid extensions
    $ python setup.py install --user # install



.. _distribution-error: 

I get an ``ImportError`` or ``DistributionError`` when using :data:`plastid`
............................................................................

If you get an error like the following::

 .. code-block:: shell

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


One or more dependencies (in this example, `scipy`_ is not installed).
Please see :ref:`this workaround <faq-install-numpy-first>`.


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


 .. toctree::
    :maxdepth: 2
