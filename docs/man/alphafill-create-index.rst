alphafill-create-index
======================

Synopsis
--------

**alphafill-create-index** [*options*]

Description
-----------

:program:`alphafill-create-index` is the program that creates a new *FastA* file for the sequences found in the PDB or PDB-REDO files.

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

.. option:: --pdb-dir=dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: --pdb-fasta=filename
	
	The FastA file containing the PDB sequences

.. option:: --threads=nr-of-threads, -t nr-of-threads
	
	Number of threads to use, zero means all available cores.

	Default is to use as many cores as the system has.

See also
--------

:manpage:`alphafill-process`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-server`, :manpage:`alphafill.conf`
