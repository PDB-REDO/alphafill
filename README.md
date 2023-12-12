[![github CI](https://github.com/pdb-redo/alphafill/actions/workflows/cmake-multi-platform.yml/badge.svg)](https://github.com/pdb-redo/alphafill/actions)
[![GitHub License](https://img.shields.io/github/license/pdb-redo/alphafill)](https://github.com/pdb-redo/alphafill/LICENSE)

# AlphaFill

**AlphaFill** is an algorithm based on sequence and structure similarity that “transplants”
missing compounds to the [AlphaFold models](https://alphafold.ebi.ac.uk/). By adding the molecular context to the protein structures, the
models can be more easily appreciated in terms of function and structure integrity.

## Building

In order to build alphafill, you need to have a modern C++ compiler (c++17), a recent version of [cmake](https://cmake.org/) and the following libraries installed:

- [Libzeep](https://github.com/mhekkel/libzeep) version 6.0.11 or higher
- [libcif++](https://github.com/PDB-REDO/libcifpp) version 6.0.0 or higher
- [libmcfp](https://github.com/mhekkel/libmcfp)
- libpq, the PostgreSQL library
- [libpqxx](http://www.pqxx.org/) version 7.2 or higher

And then you also need [yarn](https://yarnpkg.com/) to package the data for the web interface. Additionally you need [mrc](https://github.com/mhekkel/mrc) to package all the runtime data into resources in the final excutable.

Once all the requirements are met, building is as simple as:

```console
git clone https://github.com/PDB-REDO/alphafill
cd alphafill
yarn  # will fetch all node modules
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
cmake --install .
```

The building method is not OS-specific and has been tested on Ubuntu LTS 20:04. Typical compilation time is in the rorder of minutes.

## Configuration

### alphafill

For alphafill you need a copy of either PDB-REDO or PDB in mmCIF format. You also need a FastA formatted file for all sequences in this databank and an af-ligands.cif file. Since version 2.1 of alphafill This FastA file needs a custom format and can be created using alphafill itself:

```bash
alphafill create-index --pdb-dir=${PDB_DIR} --pdb-fasta=${PDB_FASTA}
```

All basic options can be stored in a file called alphafill.conf which can be located in either the current working directory or in the directory `.config` in your home directory.

Running alpafill is then as easy as:

```bash
alphafill process /srv/data/afdb/cif/AF-XXX.cif.gz /srv/data/af-filled/AF-XXX.cif.gz
```

Typical running time is less than 2 minutes but varies depending on the number of transplants.

### web interface

Before running the web application, you need to create a PostgreSQL database first. The owner and name can then be recorded in a af-filledd.conf file similar to the config file of alphafill. The location of your filled structures should be recorded as well as the db-dir option.

The database can be filled with the `alphafill rebuild-db` command. If you built the code with resources using [mrc](https://github.com/mhekkel/mrc) this will take care of setting up the tables as well, otherwise you have to create the tables running the db-schema.sql file.

After this setting up, you can start a web server using `alphafill server start`. Use `alphafill server status` to find the status of the server and `alphafill server stop` to stop it again. In this case the alphafill server runs as a daemon and log files will be written to /var/log/alphafill/. You can also start with the extra --no-daemon option and then the server will run in the foreground.
