alphafill-rebuild-db
====================

Synopsis
--------

**alphafill-rebuild-db** [*options*]

Description
-----------

:program:`alphafill-rebuild-db` is the program that creates a PostgreSQL database for the files processed by :program:`alphafill-process`. This database is required for the web application.

This program will drop an already existing database and will recreate a new one and fill it with the information found. The `db-*` options are not required, the default user credentials are used if they are not specified.

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

.. option:: --db-dir=dirname
	
	Directory containing the alphafilled data

.. option:: --pdb-dir=dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: --structure-name-pattern=value
	
	Pattern for locating structure files.

	Default is ``${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.gz``

.. option:: --metadata-name-pattern=value
	
	Pattern for locating metadata files.

	Default is ``${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.json``

.. option:: --pdb-name-pattern=value
	
	Pattern for locating PDB files.

	Default is ``${pdb-dir}/${id:1:2}/${id}/${id}_final.cif``

.. option:: --db-dbname=name
	
	AlphaFill database name

.. option:: --db-user=name
	
	AlphaFill database owner

.. option:: --db-password=value
	
	AlphaFill database password

.. option:: --db-host=value
	
	AlphaFill database host

.. option:: --db-port=value
	
	AlphaFill database port

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-process`, :manpage:`alphafill-server`, :manpage:`alphafill.conf`
