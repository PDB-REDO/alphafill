alphafill-create-index
======================

Synopsis
--------

**alphafill-create-index** [*options*]

Description
-----------

:program:`alphafill-create-index` is the program that creates a new *FastA* file for the sequences found in the PDB or PDB-REDO files.

Options
-------

.. program:: alphafill-create-index

.. option:: pdb-dir dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: pdb-fasta filename
	
	The FastA file containing the PDB sequences

.. option:: threads,t nr-of-threads
	
	Number of threads to use, zero means all available cores.

	Default is to use as many cores as the system has.

See also
--------

:manpage:`alphafill-process`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-server`
