Contributing
============

We welcome contributions! But we are new at this, so please be patient. Right
now, we (try to) follow these conventions:

Workflow
--------

.. TODO : update all repository links

Start with a ticket
...................
Before submitting a patch, please open a ticket. This will allow public
discussion of the best way to solve an issue or design a feature, and allows us
to keep logical records of everything.


Create a fork
.............
Rather than push changes to our main repository, create a fork, and, if appropriate,
a topic branch. Build & test your changes there, merge into the master branch of
your fork, then submit a pull request.


Testing
.......
We advocate test-driven development. Feature additions will not be accepted without
companion tests, and, where appropriate, test datasets. Before submitting a change,
please:

 #. Write new tests. Consider using the existing test dataset.

 #. Run your new tests, make sure they pass

 #. Run the existing tests.
       
      - If they pass, great!

      - If they fail because of an intentional design change in an object's behavior,
        make corresponding changes to the test dataset and discuss this heavily
        in the update- this will probably require a version bump.

      - If they fail otherwise, fix your submission so the old tests pass.


Finally, before submission
..........................
 #. Document everything heavily, updating module & object docstrings where
    appropriate, as well as any :doc:`/examples` that might be affected
    by the change or addition.

 #. If you use technical terms, please check if a synonym of your term is already defined
    in the :doc:`glossary <glossary>`, and then use that. If no synonym is present, please
    add the :doc:`glossary <glossary>`, and refer to it using the ``:term:`` directive.


Document formatting
-------------------
Code should be formatted as described in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_,
especially noting that we use four spaces for indentation.


Docstring & module documenation
...............................
Docstrings & documentation should be formatted in `reStructuredText`_ using
the most human-readable raw form possible. This means:

  - Format docstrings for `numpydoc`_,  as described in the
    `numpy documentation guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_

  - Follow `PEP257 <https://www.python.org/dev/peps/pep-0257>`_. for docstring style

  - Refer to classes defined in :data:`yeti` in docstrings using the shortcut `\|ClassName\|`
    rather than as `:py:class:\`package.module.ClassName\``, because the shortcuts are 
    easier to read.

  - Similarly, refer to command-line scripts as `\|modname\|`, where `modname`
    is the name of the module, excluding the extension. e.g. `\|metagene\|` for
    `:mod:\`yeti.bin.metagene\``.

  - Decorate any function that you don't want to be put in the online
    documentation with the :func:`~yeti.util.services.decorators.skipdoc`
    decorator

  - Module docstrings for command-line scripts should have a section called
    *Output files*, and, where appropriate, examples of how to call the script. 



