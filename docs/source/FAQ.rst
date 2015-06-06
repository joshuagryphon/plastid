Frequently-asked questions
==========================

.. _install_run_faq:

Installation and runtime errors
-------------------------------

.. _install_fails_virtualenv:

Install fails in a Python `Virtualenv`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   This is a due to a `known bug <https://github.com/numpy/numpy/issues/2434>`_ 
   with `NumPy`_, `SciPy`_, and `matplotlib`_ installation under setuptools. There is,
   however, a workaround. Enter your `Virtualenv`_ and manually install the SciPy
   stack via `Pip`_::

       (venv)$ pip install numpy scipy matplotlib pandas
       (venv)$ easy_install yeti.egg --reinstall


   Then check your installations by running unit tests::

       (venv)$ python -c "import numpy; numpy.test()"
       (venv)$ python -c "import scipy; scipy.test()"


   And make sure other packages are just importable::

       (venv)$ python -c "import matplotlib"
       (venv)$ python -c "import pysam"
   
   
   Then repeat the installation::
   
       (venv)$ pip install yeti


.. _distribution-error: 

I get ``ImportError`` and/or ``DistributionError`` s when using ``yeti``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


    This is due to a known problem with `NumPy`_ and `SciPy`_ dependency parsing
    during installation. See the workaround in :ref:`install_fails_virtualenv`.



    
.. _analysis_faq:
 
Analysis
--------

How do I get started?
^^^^^^^^^^^^^^^^^^^^^
pass


.. _analysis_fractional_counts:

Why do some scripts report fractional count numbers?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
pass


How do I prepare output for `DESeq`_?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
pass


