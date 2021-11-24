AlphaFill
=========

**AlphaFill** is an algorithm based on sequence and structure similarity that “transplants”
missing compounds to the [AlphaFold models](https://alphafold.ebi.ac.uk/). By adding the molecular context to the protein structures, the
models can be more easily appreciated in terms of function and structure integrity.

Building
--------
In order to build alphafill, you need to have a modern C++ compiler (c++17), a recent version of [cmake](https://cmake.org/) and the following libraries installed:

- [Libzeep](https://github.com/mhekkel/libzeep) version 5.1.5 or higher
- [libcif++](https://github.com/PDB-REDO/libcifpp) version 3.0.0 or higher (currently that's the develop branch)
- libpq, the PostgreSQL library
- [libpqxx](http://www.pqxx.org/) version 7.2 or higher

And then you also need [yarn](https://yarnpkg.com/) to package the data for the web interface. Additionally you need [mrs](https://github.com/mhekkel/mrs) to package all the runtime data into resources in the final excutable.

Once all the requirements are met, building is as simple as:

```
git clone https://github.com/PDB-REDO/alphafill
cd alphafill
yarn		# will fetch all node modules
mkdir build
cd build
cmake ..
cmake --build .
cmake --install .
```

