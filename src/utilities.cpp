/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 Maarten L. Hekkelman, NKI-AVL
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <filesystem>
#include <fstream>
#include <iostream>

#include <cif++.hpp>
#include <cif++/compound.hpp>

#include "revision.hpp"
#include "utilities.hpp"

namespace fs = std::filesystem;
namespace po = boost::program_options;

// --------------------------------------------------------------------

std::unique_ptr<file_locator> file_locator::s_instance;

void file_locator::init(boost::program_options::variables_map &vm)
{
	fs::path dbDir = vm["db-dir"].as<std::string>();
	if (not fs::is_directory(dbDir))
		throw std::runtime_error("AlphfaFill data directory does not exist");

	fs::path pdbDir = vm["pdb-dir"].as<std::string>();
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	s_instance.reset(new file_locator(dbDir, pdbDir, vm["structure-name-pattern"].as<std::string>(),
		vm["pdb-name-pattern"].as<std::string>(), vm["metadata-name-pattern"].as<std::string>()));
}

std::vector<std::filesystem::path> file_locator::get_all_structure_files(const std::string &id)
{
	std::vector<fs::path> result;

	int i = 1;
	for (;;)
	{
		fs::path chunk = s_instance->get_structure_file(id, i);
		if (not fs::exists(chunk))
			break;
		
		result.emplace_back(std::move(chunk));
		++i;
	}
	
	return result;
}

// --------------------------------------------------------------------

po::variables_map load_options(int argc, char *const argv[],
	po::options_description &visible_options,
	po::options_description &hidden_options, po::positional_options_description &positional_options,
	const char *config_file_name)
{
	visible_options.add_options()
		("af-dir", po::value<std::string>(), "Directory containing the alphafold data")
		("db-dir", po::value<std::string>(), "Directory containing the af-filled data")
		("pdb-dir", po::value<std::string>(), "Directory containing the mmCIF files for the PDB")

		("ligands", po::value<std::string>()->default_value("af-ligands.cif"), "File in CIF format describing the ligands and their modifications")

		("compounds", po::value<std::string>(), "Location of the components.cif file from CCD")
		("components", po::value<std::string>(), "Location of the components.cif file from CCD, alias")
		("extra-compounds", po::value<std::string>(), "File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file")
		("mmcif-dictionary", po::value<std::string>(), "Path to the mmcif_pdbx.dic file to use instead of default")

		("structure-name-pattern", po::value<std::string>(), "Pattern for locating structure files")
		("metadata-name-pattern", po::value<std::string>(), "Pattern for locating metadata files")
		("pdb-name-pattern", po::value<std::string>(), "Pattern for locating PDB files")

		("config", po::value<std::string>(), "Config file")

		("help,h", "Display help message")
		("version", "Print version")
		("verbose,v", "Verbose output")
		("quiet", "Do not produce warnings");

	hidden_options.add_options()
		("debug,d", po::value<int>(), "Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv)
				  .options(cmdline_options)
				  .positional(positional_options)
				  .run(),
		vm);

	fs::path configFile = config_file_name;
	if (vm.count("config"))
		configFile = vm["config"].as<std::string>();

	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / configFile;

	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options, true), vm);
	}

	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		write_version_string(std::cout, vm.count("verbose"));
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cout << visible_options << std::endl;
		exit(0);
	}

	if (vm.count("quiet"))
		cif::VERBOSE = -1;

	if (vm.count("verbose"))
		cif::VERBOSE = 1;

	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	if (vm.count("db-dir") == 0)
	{
		std::cout << "AlphaFill data directory not specified" << std::endl;
		exit(1);
	}

	if (vm.count("pdb-dir") == 0)
	{
		std::cout << "PDB directory not specified" << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------
	// Load extra CCD definitions, if any

	if (vm.count("compounds"))
		cif::add_file_resource("components.cif", vm["compounds"].as<std::string>());
	else if (vm.count("components"))
		cif::add_file_resource("components.cif", vm["components"].as<std::string>());

	if (vm.count("extra-compounds"))
		cif::compound_factory::instance().push_dictionary(vm["extra-compounds"].as<std::string>());

	// And perhaps a private mmcif_pdbx dictionary

	if (vm.count("mmcif-dictionary"))
		cif::add_file_resource("mmcif_pdbx_v50.dic", vm["mmcif-dictionary"].as<std::string>());

	return vm;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

int main(int argc, char *const argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::add_data_directory(DATA_DIR);
#endif
		result = a_main(argc, argv);
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
