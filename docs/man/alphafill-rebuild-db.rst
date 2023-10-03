alphafill-rebuild-db
=================

Synopsis
--------

**alphafill-rebuild-db** [*options*]

Description
-----------

:program:`alphafill-rebuild-db` is the program that creates a PostgreSQL database for the files processed by :program:`alphafill-process`. This database is required for the web application.

This program will drop an already existing database and will recreate a new one and fill it with the information found. The `db-*` options are not required, the default user credentials are used if they are not specified.

Options
-------

.. program:: alphafill-rebuild-db

.. option:: db-dir dirname
	
	Directory containing the alphafilled data

.. option:: pdb-dir dirname
	
	Directory containing the mmCIF files for the PDB

.. option:: structure-name-pattern
	
	Pattern for locating structure files.

	Default is `${db-dir}/${id:0:2}/AF-${id}-F${chunk}-model_v${version}.cif.gz`

.. option:: metadata-name-pattern
	
	Pattern for locating metadata files.

	Default is `${db-dir}/${id:0:2}/AF-${id}-F${chunk}-model_v${version}.cif.json`

.. option:: pdb-name-pattern
	
	Pattern for locating PDB files.

	Default is `${pdb-dir}/${id:1:2}/${id}/${id}_final.cif`

.. option:: pae-name-pattern
	
	Pattern for location cached PAE scores.

	Default is `${db-dir}/${id:0:2}/AF-${id}-F${chunk}-model_v${version}.pae.json`

.. option:: db-dbname
	
	AlphaFill database name

.. option:: db-user
	
	AlphaFill database owner

.. option:: db-password
	
	AlphaFill database password

.. option:: db-host
	
	AlphaFill database host

.. option:: db-port
	
	AlphaFill database port

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-process`, :manpage:`alphafill-server`
