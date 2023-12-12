alphafill-server
=================

Synopsis
--------

**alphafill-server** [*options*] <*command*>

Description
-----------

:program:`alphafill-server` is the web application server for alphafill. This is a single executable containing everything required to serve AlphaFill as a web service.

Before running the server, you should have created some AlphaFill files using :program:`alphafill-process` and created a PostgreSQL database with :program:`alphafill-rebuild-db`.

The **command** argument must be one of:

start
	Starts a new server process. The program is usually run as a daemon process but for debugging purposes you can run it in the foreground using the :option:`--no-daemon,-F` flag.

	If run in the background, a separate process is forked and will drop privileges while the main process listening to the :option:`--port` is run using the credentials of the user starting the process. Each connection is passed to child process for handling.

	When running in the background, two log files are written to ``/var/log/alphafill``, one called *access.log* and the other *error.log*. And a file is written to ``/var/run/alphafill`` containing the process ID of the daemon process.

status
	Report the status of a currently running server instance. Result can be either *running* or *stopped*.

stop
	Stop the currently running server process.

reload
	Stops the child process and closes all log files. Then restarts a new child process.

The options for the server are as follows.

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

.. option:: --no-daemon,-F
	
	Do not fork a background process but run the web server in the foreground.
	
.. option:: --address=value
	
	Address to listen to.
	
	Default value is *127.0.0.1* (i.e. localhost)
	
.. option:: --port=value
	
	Port to listen to.

	Default value is *10342*
	
.. option:: --user=name
	
	User to run as.

	Default value is *www-data*
	
.. option:: --context=value
	
	Reverse proxy context.

	When the server is supposed to be accessible from the outside, you'd best put a reverse proxy server before it since HTTPS is not supported. If you do so, the external address can be provided in this option to generate correct links in the web pages.
	
.. option:: --db-link-template=value
	
	Template for links to *PDB* or *PDB-REDO* entries. Result pages contain PDB-IDs that have a link. To make them point to something outside the scope of alphafill, you can provide a link template in this option. The *variable* ``${id}`` will be replaced with the PDB-ID referenced.

.. option:: --db-dbname=name
	
	The name of the AlphaFill PostgreSQL database.

.. option:: --db-user=name
	
	The owner of the AlphaFill PostgreSQL database.

.. option:: --db-password=value
	
	The password of the AlphaFill PostgreSQL database.

.. option:: --db-host=value
	
	The host of the AlphaFill PostgreSQL database.

.. option:: --db-port=value
	
	The port of the AlphaFill PostgreSQL database.

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

.. option:: --blast-no-filter
	
	By default blast will use a low complexity filter. Use this option to turn that off.	

.. option:: --blast-no-gapped
	
	By default blast performs gapped alignment. Use this option to turn that off.

.. option:: --blast-gap-open=value
	
	Blast penalty for gap open.

	Default value is 11.

.. option:: --blast-gap-extend=value
	
	Blast penalty for gap extend.

	Default value is 1.

.. option:: --threads=value, -t value
	
	Number of threads to use, zero means all available cores.

	Default is to use as many cores as the system has.

.. option:: --structure-name-pattern=value
	
	Template used for locating structure files.

	Default value is ``${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.gz``

.. option:: --metadata-name-pattern=value
	
	Template used for locating metadata files
	
	Default value is ``${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.json``

.. option:: --pdb-name-pattern=value
	
	Template used for locating PDB files

	Default value is ``${pdb-dir}/${id:1:2}/${id}/${id}_final.cif``
	
.. option:: --alphafold-3d-beacon=value
	
	The URL of the 3d-beacons service for alphafold

	Default value is ``https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/${id}.json?provider=alphafold``

.. option:: --custom-dir=dirname
	
	Directory for custom built entries. These are files uploaded by the user of the web service.

	Default value is ``/tmp/alphafill``

.. option:: --yasara=filename
	
	Location of the yasara executable, needed for optimising.

	Default value is ``/opt/yasara/yasara``

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-process`, :manpage:`alphafill.conf`
