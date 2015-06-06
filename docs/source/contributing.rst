Contributing
============

We haven\'t thought much about this yet, but at the moment we're adhering to several conventions:

1. We follow `Vincent Driessen's git branching model <http://nvie.com/posts/a-successful-git-branching-model/>`_.
   `git-flow <https://github.com/nvie/gitflow/wiki/Installation>`_
   provides tools useful for this.
       
2. Code should be formatted as described
   in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_,
   especially noting that we use four spaces for indentation.

3. Docstrings should be formatted as described using
   the `numpydoc plugin <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
   for `Sphinx <http://sphinx-doc.org/>`_,
   with guidance from `PEP257 <https://www.python.org/dev/peps/pep-0257>`_.
   This means that docstrings are formatted in
   `reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_.
   Beyond that, please reference classes in docstrings using the shortcut `|ClassName|`
   rather than `:py:class:\`package.module.ClassName\``, because we perform these substitutions
   automatically and generate cross-references in the doc build process. It also makes
   docstrings more readable to humans perusing the code.

