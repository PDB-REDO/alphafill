alphafill-process
=================

Synopsis
--------

**alphafill-process** [*options*] <*inputfile*> [<*outputfile*>]

Description
-----------

:program:`alphafill-process` is the program to create a new AlphaFill file based on an input model.

.. option:: --help

	Display the options allowed for this program.

.. option:: --version

	Display the version of this program.

.. option:: --verbose

	Use a more verbose output, printing status and progress information.

.. option:: --quiet

	Do not print any status or progress information.

.. option:: --config configfile

	Use the file *configfile* to collection options. The default is to look for a file called *alphafill.conf* in the current directory and then in the directory */etc*. Use this option to override this and specify your own configuration file.

.. option:: --pae-file=filename
	
	Specify a specific file containing PAE information, default is to use a filename based on inputfile

.. option:: --db-dir=dirname
	
	Directory containing the alphafilled data

.. option:: --pdb-dir=dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: --pdb-fasta=filename
	
	The FastA file containing the PDB sequences

.. option:: --ligands=filename
	
	File in CIF format describing the ligands and their modifications.
	
	The default file is af-ligands.cif	

.. option:: --max-ligand-to-backbone-distance=value
	
	The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure.
	
	Default value is 6.	

.. option:: --min-hsp-identity=value
	
	The minimal identity for a high scoring pair (note, value between 0 and 1).
	
	Default value is 0.25.

.. option:: --min-alignment-length=value
	
	The minimal length of an alignment.

	Default value is 85.	

.. option:: --min-separation-distance=value
	
	The centroids of two identical ligands should be at least this far apart to count as separate occurrences.

	Default value is 3.5.

.. option:: --clash-distance-cutoff=value
	
	The max distance between polymer atoms and ligand atoms used in calculating clash scores.

	Default value is 4.

.. option:: --blast-report-limit=value
	
	Number of blast hits to use.

	Default value is 250.	

.. option:: --blast-matrix=value
	
	Blast matrix to use.

	Default matrix is *BLOSUM62*.

.. option:: --blast-word-size=value
	
	Blast word size.

	Default value is 3.

.. option:: --blast-expect=value
	
	Blast expect cut off.

	Default value is 10.

.. option:: --blast-no-filter=value
	
	By default blast will use a low complexity filter. Use this option to turn that off.	

.. option:: --blast-no-gapped
	
	By default blast performs gapped alignment. Use this option to turn that off.

.. option:: --blast-gap-open
	
	Blast penalty for gap open.

	Default value is 11.

.. option:: --blast-gap-extend=value
	
	Blast penalty for gap extend.

	Default value is 1.

.. option:: --threads=value, -t value
	
	Number of threads to use, zero means all available cores.

	Default is to use as many cores as the system has.

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-server`, :manpage:`alphafill.conf`
