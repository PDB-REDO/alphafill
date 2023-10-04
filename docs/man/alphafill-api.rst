alphafill-api
=============

Description
-----------

The :program:`alphafill-server` application, when running, can be accessed using REST calls.

The url's below are specifying https://alphafill.eu as server. 

GET \https://alphafill.eu/v1/aff/{id}
	Retrieve the AlphaFill structure for ID ``{id}``.

GET \https://alphafill.eu/v1/aff/{id}/json
	Retrieve the accompanying JSON file for ID ``{id}``.

GET \https://alphafill.eu/v1/aff/{id}/status
	Retrieve the status for ID ``{id}``.

GET \https://alphafill.eu/v1/aff/3d-beacon/{id}
	Retrieve the 3d-beacon information for ID ``{id}``.

GET \https://alphafill.eu/v1/aff/{id}/stripped/{asymlist}
	Retrieve a stripped down structure file for ID ``{id}`` where only the asyms in ``{asymlist}`` are included.

GET \https://alphafill.eu/v1/aff/{id}/stripped/{asymlist}/{identity}
	Retrieve a stripped down structure file for ID ``{id}`` where only the asyms in ``{asymlist}`` are included using an identity cut-off of ``{identity}``.

GET \https://alphafill.eu/v1/aff/{id}/optimized/{asymlist}
	Retrieve an optimized version of the structure file for ID ``{id}`` where the asyms in ``{asymlist}`` and their environment are optimized using yasara.

GET \https://alphafill.eu/v1/aff/{id}/optimized-with-stats/{asymlist}
	Retrieve an optimized version of the structure file for ID ``{id}`` where the asyms in ``{asymlist}`` and their environment are optimized using yasara. This version also includes statistical data.

POST \https://alphafill.eu/v1/aff
	Upload a structure for processing by *AlphaFill*. The body should contain a parameter called ``structure`` that has as value the AlphaFold file you want to process. An optional parameter ``pae`` may contain the *PAE* score matrix data for this file.

	The result of this call will be a JSON object containing a field ``id`` with the assigned ID, or a field name ``error`` describing what went wrong.

See also
--------

:manpage:`alphafill-create-index`, :manpage:`alphafill-process`, :manpage:`alphafill-rebuild-db`, :manpage:`alphafill-process`
