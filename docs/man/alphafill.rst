alphafill
=========

Synopsis
--------

The alphafill executable takes as first argument a program to run. The available programs are:

**alphafill create-index**
	Tool to create a *FastA* file for the PDB or PDB-REDO structure files

**alphafill process**
	Tool to process an AlphaFold structure

**alphafill rebuild-db**
	Tool to create a database with processed structures, used by the web application

**alphafill server**
	Tool to serve the processed entries as a web application

Description
-----------

*AlphaFill* is an algorithm based on sequence and structure similarity that “transplants” missing ligands, cofactors and (metal) ions to the AlphaFold models. By adding the molecular context to these protein structures, the models can be more easily appreciated in terms of function and structural integrity. Consequently, the *AlphaFill* models can be helpful in designing downstream wet-lab experiments and/or computational studies.

The PDB files used by *AlphaFill* are normally `PDB-REDO <https://pdb-redo.eu>`_ files. These files should be located in a directory specified by the :option:`--db-dir` option. The default for the :option:`--pdb-name-pattern` option assumes you are using these *PDB-REDO* files. If you opt to use another source of PDB files, you will have to change the value for this option.

E.g. for regular PDB files, the value should be: ``${pdb-dir}/${id:1:2}/${id}.cif.gz``


All programs accept at least the following options

.. option:: --help

	Display the options allowed for this program.

.. option:: --version

	Display the version of this program.

.. option:: --verbose

	Use a more verbose output, printing status and progress information.

.. option:: --quiet

	Do not print any status or progress information.

.. option:: --config=configfile

	Use the file *configfile* to collection options. The default is to look for a file called *alphafill.conf* in the current directory and then in the directory */etc*. Use this option to override this and specify your own configuration file.

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-process`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-process`
