Extending plastid
=================

Plastid defines the following :term:`entry points` to enable users to write
plug-in functions that can be accessed from the command line:


    ===========================    =================================================== 
     **Entry point**                **Used for**
    ---------------------------    --------------------------------------------------- 

    ``plastid.mapping_rules``      Adding new :term:`mapping rules <mapping rule>`
                                   to :data:`plastid`'s command-line scripts.

    ``plastid.mapping_options``    Adding command-line options used by new
                                   :term:`mapping rules <mapping rule>`.

    ===========================    =================================================== 


For :data:`plastid` to discover your plug-in, your plug-in must be registered
with the system. Registration requires some packaging, which isn't too painful.
Packaging is discussed in the following sections:

  - :ref:`entrypoints-folder`
  - :ref:`entrypoints-write-options`

      - :ref:`entrypoints-mapping-functions`
      - :ref:`entrypoints-parameters`

  - :ref:`entrypoints-setup-py`

  - :ref:`entrypoints-install`


 .. _entrypoints-folder:


Setting up the folder
---------------------
First, create a folder with structure similar to the following:

 .. code-block:: none

    my_project/
        setup.py
        my_project/
            __init__.py
            map_rules.py


Adjust the filenames to suit your project.


 .. _entrypoints-write-options:



Writing command-line options for mapping rules
----------------------------------------------

 .. _entrypoints-mapping-functions:



Mapping rules
.............

We assume you have written mapping rules as described in 
:ref:`mapping-rules-roll-your-own`. :data:`plastid` needs some metadata
to use them. This is specified in a dictionary that defines at least
`bamfunc` or `bowtiefunc`. All of the remaining keys are optional:

    ====================  =================  =============================================
     **Key**              **Value type**     **Value**
    --------------------  -----------------  ---------------------------------------------

    `name`                str                Overrides the command-line name
                                             of mapping rule defined in ``setup.py``.
                                             I.e. - the flag ``--name`` will be the
                                             command-line argument that invokes the rule

                                             Must not contain spaces, dashes, or special
                                             characters. Underscores are o

    `bamfunc`             Function           Mapping function for alignments
                                             in `BAM`_ format

    `bowtiefunc`          Function           Mapping function for alignments
                                             in `bowtie`_ format

    `help`                str                Command-line help for the 
                                             mapping function. Should describe
                                             what the function does, and 
                                             which command-line arguments
                                             affect its behavior (e.g. 
                                             ``--offset``, ``--nibble`` or
                                             something added in 
                                             :ref:`entrypoints-parameters`)

    ====================  =================  =============================================


If `bowtiefunc` or `bamfunc` are unspecified or set to `None`, 
:data:`plastid` will assume the mapping function is not implemented 
for the corresponding type. Typically, users would only write
a function for mapping `BAM`_ files.


We'll suppose that all of our functions are specified in ``my_project/map_rules.py``
as described :ref:`above <entrypoints-folder>`. The contents of ``map_rules.py``
might then look something like this:

 .. code-block:: python

    #!/usr/bin/env python

    def rule1_for_bowtie_files(alignment,args=None):
        # calculate position(s) where a single aliignment maps
        # and the value to place at each position
        #
        # the parsed command-line arguments will be passed
        # as an argparse.Namespace object
        ...

        return position_value_tuples

    def rule1_for_BAM_files(alignments,segment,args=args):
        # calculate positions where a list of alignments map,
        # and a vector of values at each position
        #
        # again, args is an argparse.Namespace object
        # from the command-line args
        ...

        return reads_out, count_array

    def rule2_for_BAM_files_only(alignments,segment,args=args):
        # calculate positions where a list of alignments map,
        # and a vector of values at each position
        ...

        # do something with a command-line argument
        my_option = args.new_option
        if my_option == "":
            pass

        return reads_out, count_array


    rule1_info = {
        "name"       : 'rule1',
        "bamfunc"    : rule1_for_BAM_files,
        "bowtiefunc" : rule1_for_bowtie_files,
        "help"       : "Some help text for rule 1."
    }


    rule2_info = {
        "name"       : 'rule2',
        "bamfunc"    : rule2_for_BAM_files_only,
        "help"       : "Some help text. Rule 2's behavior is modified by the option `--new_option`"
    }


`rule1` is defined for both `BAM`_ and `bowtie`_ files. `rule2` is defined
only for `BAM`_ files, and it uses the command-line option ``--new_option``,
which we define below in :ref:`entrypoints-parameters`.


 .. _entrypoints-parameters:

Additional parameters for mapping rules
.......................................

Additional command-line parameters are also specified as dictionaries.
In these, the keys and values can be any valid parameters for
:meth:`argparse.ArgumentParser.add_argument`. Each dictionary should
additionally define a key called `name`, whose value will be used as
the name of the command-line argument. For example, we might add
the following lines to ``my_project/map_rules.py``:

 .. code-block:: python

    param1 = { 
        "name"  : "new_option",
        "type"  : int,
        "nargs" : 2,
        "help"  : "Some help text for --new_option",
        "metavar" : "N",
    }


That's it!



 .. _entrypoints-setup-py:

Writing ``setup.py``
--------------------

Having written the mapping functions and made dictionaries describing them,
we need to write package metadata so that :data:`plastid` can find the new
functions. All of this information goes into ``setup.py``. 

``setup.py`` should everything needed to set up and install your package.
For more information see the documentation for :mod:`setuptools` and / or
:mod:`distutils`. ``setup.py`` should minimally contain the following:

 .. code-block:: python

    #/usr/bin/env python
    from setuptools import setup, find_packages


    # list all the rules we want to include
    # syntax is: 
    #
    #    rule_name = path.to.rule:rule_info_dictionary"
    #
    #
    rules = [
        "rule1 = my_project.rules:rule1_info",
        "rule2 = my_project.rules:rule2_info",
    ]

    # list any extra arguments we want to include
    # syntax is: 
    #
    #    argument_name = path.to.rule:arg_info_dictionary"
    #
    #
    rule_options = [
        "new_option = my_project.rules:param1",
    ]


    setup(
        # root level name of package
        name = "my_project",

        # tell setup() that `rules` and `rule_options` specify mapping
        # ruls and arguments for plastid:
        entry_points = { 
            "plastid.mapping_rules"   : rules,
            "plastid.mapping_options" : rule_options,
        },

        setup_requires = ['plastid>=0.4.4'],
        packages = find_packages(),

        # plus any other arguments (e.g. package author, description)
        # to ``setup``. 

    )

That's the last piece.

 
  .. _entrypoints-install:

Installing the new mapping rules
--------------------------------

Installation is the final step. Enter the folder containing ``setup.py``. 
Then, to install your new mapping rules, type:

 .. code-block:: shell

    $ python setup.py install [--user]

 .. 
 
 
Or, if you plan to keep developing your :term:`mapping rules <mapping rule>`,
and want :data:`plastid` to be aware of these changes instantly:

 .. code-block:: shell

    $ python setup.py develop --user


To test your installation, check command-line help from a script that uses
mapping rules (e.g. ``make_wiggle``):

 .. code-block:: shell

    $ make-wiggle --help

If the installation proceeded correctly you should see something like this:

 .. code-block:: none

    # rest of command line help above


    alignment mapping options (BAM & bowtie files only):
      For BAM or bowtie files, one of the mutually exclusive read mapping choices
      is required:

      --fiveprime_variable  Map read alignment to a variable offset from 5'
                            position of read, with offset determined by read
                            length. Requires `--offset` below
      --fiveprime           Map read alignment to 5' position.
      --threeprime          Map read alignment to 3' position
      --center              Subtract N positions from each end of read, and add
                            1/(length-N), to each remaining position, where N is
                            specified by `--nibble`
      --rule2               Some help text. Rule 2's behavior is modified by the
                            option `--new_option`
      --rule1               Some help text for rule 1.

      
      The remaining arguments are optional and affect the behavior of specific
      mapping rules:

      --offset OFFSET       For `--fiveprime` or `--threeprime`, provide an
                            integer representing the offset into the read,
                            starting from either the 5' or 3' end, at which data
                            should be plotted. For `--fiveprime_variable`, provide
                            the filename of a two-column tab-delimited text file,
                            in which first column represents read length or the
                            special keyword `'default'`, and the second column
                            represents the offset from the five prime end of that
                            read length at which the read should be mapped.
      --nibble N            For use with `--center` only. nt to remove from each
                            end of read before mapping (Default: 0)
      --new_option N N      Some help text for --new_option


    # remaining command-line help below


If the new mapping rule and command-line arguments are listed, you are ready.



------------------------------------------------------------------------------

See also
--------

  - :doc:`/concepts/mapping_rules` for information on how to write
    :term:`mapping rules <mapping rule>`

  - :mod:`argparse` documentation for information on command-line arguments

  - Documentation for :mod:`setuptools` and :mod:`distutils` for more information
    on packaging

