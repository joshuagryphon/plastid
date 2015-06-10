Contributing
============

We welcome contributions! But we are new at this, so please be patient. Right
now, we follow these conventions:

1.  *Testing:* We require test-driven development. Feature additions will not
    be accepted without companion tests, and, where appropriate, test datasets.
    
2.  *Branching:* We follow `Vincent Driessen's git branching model <http://nvie.com/posts/a-successful-git-branching-model/>`_.
    `git-flow <https://github.com/nvie/gitflow/wiki/Installation>`_
    provides tools useful for this.
       
3.  *Code formatting:* Code should be formatted as described
    in `PEP8 <https://www.python.org/dev/peps/pep-0008>`_, especially noting
    that we use four spaces for indentation.

4.  *Docstring and documentation:* These should be formatted for
    `Sphinx <http://sphinx-doc.org/>`_, as described in the
    `numpy documentation guide <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
    with guidance from `PEP257 <https://www.python.org/dev/peps/pep-0257>`_.
    This means that docstrings are formatted in
    `reStructuredText <http://docutils.sourceforge.net/docs/user/rst/quickref.html>`_.

   Beyond that:
   
      1. Refer to classes in docstrings using the shortcut `\|ClassName\|`
         rather than as `:py:class:\`package.module.ClassName\``, because the
         former is easier to read and handled in post-processing
      
      2. If you use technical terms, please check if a synonym of your term
         is already defined in the :doc:`glossary <glossary>`, and then use that.
         If no synonym is present, please add yours to the
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