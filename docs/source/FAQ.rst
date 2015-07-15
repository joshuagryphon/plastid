Frequently-asked questions
==========================


-------------------------------------------------------------------------------

.. _faq-run:

Installation and runtime
------------------------

.. _faq-install-fails-virtualenv:

Install fails in a Python `Virtualenv`_
.......................................

This is a due to a `known bug <https://github.com/numpy/numpy/issues/2434>`_ 
with `NumPy`_, `SciPy`_, and `matplotlib`_ installation under setuptools. There is,
however, a workaround. Enter your `Virtualenv`_ and manually install the SciPy
stack via `Pip`_:

.. code-block:: shell

    (venv)$ pip install numpy scipy matplotlib pandas
    (venv)$ easy_install yeti.egg --reinstall


Then check your installations by running unit tests:

.. code-block:: shell

    (venv)$ python -c "import numpy; numpy.test()"
    (venv)$ python -c "import scipy; scipy.test()"


And make sure other packages are just importable:

.. code-block:: shell

    (venv)$ python -c "import matplotlib"
    (venv)$ python -c "import pysam"


Then repeat the installation:

.. code-block:: shell

    (venv)$ pip install yeti


.. _distribution-error: 

I get an ``ImportError`` or ``DistributionError`` when using :data:`yeti`
.........................................................................

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


One or more dependencies (in this example, `scipy`_ is not installed). If
installing in a `virtualenv`_, please see
:ref:`this workaround <faq-install-fails-virtualenv>`.


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


.. _faq-cs-vs-counts-in-region:

What is the difference between :mod:`~yeti.bin.counts_in_region` and :mod:`~yeti.bin.cs`?
.........................................................................................
:mod:`~yeti.bin.counts_in_region` very simply counts read coverage (or any data) over
regions of interest, and reports those numbers in terms of :term:`counts` and :term:`RPKM`. It can 
optionally take a :term:`mask file`, if there are genomic positions in the regions
of interest which should be excluded from analysis. Otherwise, it makes no corrections.

:mod:`~yeti.bin.cs` is more complex, and is principally designed to make rough estimates
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

Either one is an appropriate starting place for a pipeline, depending upon your needs.


.. _faq-analysis-deseq:

How do I prepare output for `DESeq`_?
.....................................

TODO: write this


 .. toctree::
    :maxdepth: 2
