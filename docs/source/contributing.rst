Contributing
============

We welcome contributions! But we are new at this, so please be patient. Right
now, we follow these conventions:

Workflow
--------
Before submitting a patch, please open a ticket. This will allow public discussion
of the best way to solve an issue or design a feature, and allows us to keep
logical records of everything.

Testing
-------
We advocate test-driven development. Feature additions will not be accepted without
companion tests, and, where appropriate, test datasets.

Document formatting
-------------------

  - Code should be formatted as described in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_,
    especially noting that we use four spaces for indentation.

  - Docstrings & documentation should be formatted in `reStructuredText`_ using
    the most human-readable raw form possible. This means:

      - Format docstrings for `numpydoc`_,  as described in the
        `numpy documentation guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

      - Follow `PEP257 <https://www.python.org/dev/peps/pep-0257>`_. for docstring style

      - Refer to classes defined in :data:`yeti` in docstrings using the shortcut `\|ClassName\|`
        rather than as `:py:class:\`package.module.ClassName\``, because the shortcuts are 
        easier to read.
      


      2. If you use technical terms, please check if a synonym of your term
         is already defined in the :doc:`glossary <glossary>`, and then use that.
         If no synonym is present, please add the
         :doc:`glossary <glossary>`, and refer to it using the ``:term:``
         directive.
         
      3. Decorate any function that you don't want to be put in the online
         documentation with the :func:`~yeti.util.services.decorators.skipdoc`
         decorator
      
      4. Command-line scripts should have a section called *Output files*
         in their module docstrings, and, where appropriate, examples of how
         to call the script. 

    We are aware that we don't always manage to adhere to the above. We are
    in the process of standardizing.
