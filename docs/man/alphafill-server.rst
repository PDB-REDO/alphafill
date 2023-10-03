alphafill-server
=================

Synopsis
--------

**alphafill-server** [*options*] <*command*>

Description
-----------

:program:`alphafill-server` is the web application server for alphafill.

Options
-------

.. program:: alphafill-server


.. option:: pae-file filename
	
	Specify a specific file containing PAE information, default is to use a filename based on inputfile

.. option:: db-dir dirname
	
	Directory containing the alphafilled data

.. option:: pdb-dir dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: pdb-fasta filename
	
	The FastA file containing the PDB sequences

.. option:: ligands filename
	
	File in CIF format describing the ligands and their modifications.
	
	The default file is af-ligands.cif	

.. option:: max-ligand-to-backbone-distance
	
	The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure.
	
	Default value is 6.	

.. option:: min-hsp-identity
	
	The minimal identity for a high scoring pair (note, value between 0 and 1).
	
	Default value is 0.25.

.. option:: min-alignment-length
	
	The minimal length of an alignment.

	Default value is 85.	

.. option:: min-separation-distance
	
	The centroids of two identical ligands should be at least this far apart to count as separate occurrences.

	Default value is 3.5.

.. option:: clash-distance-cutoff
	
	The max distance between polymer atoms and ligand atoms used in calculating clash scores.

	Default value is 4.

.. option:: blast-report-limit
	
	Number of blast hits to use.

	Default value is 250.	

.. option:: blast-matrix
	
	Blast matrix to use.

	Default matrix is *BLOSUM62*.

.. option:: blast-word-size
	
	Blast word size.

	Default value is 3.

.. option:: blast-expect
	
	Blast expect cut off.

	Default value is 10.

.. option:: blast-no-filter
	
	By default blast will use a low complexity filter. Use this option to turn that off.	

.. option:: blast-no-gapped
	
	By default blast performs gapped alignment. Use this option to turn that off.

.. option:: blast-gap-open
	
	Blast penalty for gap open.

	Default value is 11.

.. option:: blast-gap-extend
	
	Blast penalty for gap extend.

	Default value is 1.

.. option:: threads,t nr-of-threads
	
	Number of threads to use, zero means all available cores.

	Default is to use as many cores as the system has.

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-process`
