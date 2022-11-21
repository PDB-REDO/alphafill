/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
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

#include <iostream>
#include <filesystem>
#include <thread>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include "alphafill.hpp"
#include "data-service.hpp"
#include "main.hpp"
#include "revision.hpp"
#include "server.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;

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

		auto &config = mcfp::config::instance();
		config.init(
			"usage: alphafill command [options]\n       (where command is one of 'server', 'process', 'validate' or 'prepare-pdb-list'",
			mcfp::make_option("version", "Show version number"),
			mcfp::make_option("verbose,v", "Show verbose output"),

			mcfp::make_option("help,h", "Display help message"),
			mcfp::make_option("quiet", "Do not produce warnings"),

			mcfp::make_option<std::string>("config", "alphafill.conf", "Configuration file to use"),

			mcfp::make_option<std::string>("af-dir", "Directory containing the alphafold data"),
			mcfp::make_option<std::string>("db-dir", "Directory containing the alphafilled data"),
			mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),

			mcfp::make_option<std::string>("pdb-fasta", "The FastA file containing the PDB sequences"),
			mcfp::make_option<std::string>("pdb-id-list", "Optional file containing the list of PDB ID's that have any of the transplantable ligands"),

			mcfp::make_option<std::string>("ligands", "af-ligands.cif", "File in CIF format describing the ligands and their modifications"),

			mcfp::make_option<float>("max-ligand-to-backbone-distance", 6, "The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure"),
			mcfp::make_option<float>("min-hsp-identity", 0.25, "The minimal identity for a high scoring pair (note, value between 0 and 1)"),
			mcfp::make_option<int>("min-alignment-length", 85, "The minimal length of an alignment"),
			mcfp::make_option<float>("min-separation-distance", 3.5, "The centroids of two identical ligands should be at least this far apart to count as separate occurrences"),
			mcfp::make_option<uint32_t>("blast-report-limit", 250, "Number of blast hits to use"),

			mcfp::make_option<float>("clash-distance-cutoff", 4, "The max distance between polymer atoms and ligand atoms used in calculating clash scores"),

			mcfp::make_option<std::string>("compounds", "Location of the components.cif file from CCD"),
			mcfp::make_option<std::string>("components", "Location of the components.cif file from CCD, alias"),
			mcfp::make_option<std::string>("extra-compounds", "File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file"),
			mcfp::make_option<std::string>("mmcif-dictionary", "Path to the mmcif_pdbx.dic file to use instead of default"),

			mcfp::make_option<std::string>("structure-name-pattern", "Pattern for locating structure files"),
			mcfp::make_option<std::string>("metadata-name-pattern", "Pattern for locating metadata files"),
			mcfp::make_option<std::string>("pdb-name-pattern", "Pattern for locating PDB files"),

			mcfp::make_option<int>("threads,t", std::thread::hardware_concurrency(), "Number of threads to use, zero means all available cores"),

			mcfp::make_hidden_option("validate-fasta", "Validate the FastA file (check if all sequence therein are the same as in the corresponding PDB files)"),
			mcfp::make_hidden_option("prepare-pdb-list", "Generate a list with PDB ID's that contain any of the ligands"),

			mcfp::make_hidden_option<std::string>("test-pdb-id", "Test with single PDB ID"),

			mcfp::make_option<std::string>("alphafold-3d-beacon", "The URL of the 3d-beacons service for alphafold"),
			mcfp::make_hidden_option<std::string>("test-af-id", ""),

			mcfp::make_option("no-daemon,F", "Do not fork a background process"),
			mcfp::make_option<std::string>("address", "Address to listen to"),
			mcfp::make_option<unsigned short>("port", "Port to listen to"),
			mcfp::make_option<std::string>("user", "User to run as"),
			mcfp::make_option<std::string>("context", "Reverse proxy context"),
			mcfp::make_option<std::string>("db-link-template", "Template for links to pdb(-redo) entry"),
			mcfp::make_option<std::string>("db-dbname", "AF DB name"),
			mcfp::make_option<std::string>("db-user", "AF DB owner"),
			mcfp::make_option<std::string>("db-password", "AF DB password"),
			mcfp::make_option<std::string>("db-host", "AF DB host"),
			mcfp::make_option<std::string>("db-port", "AF DB port"),

			mcfp::make_option<std::string>("custom-dir", (fs::temp_directory_path() / "alphafill").string(), "Directory for custom built entries"),

			mcfp::make_option<std::string>("yasara", "/opt/yasara/yasara", "Location of the yasara executable, needed for optimising"),

			mcfp::make_hidden_option("test", "Run test code")
			);

		// config.set_ignore_unknown(true);
		config.parse(argc, argv);

		if (config.has("version"))
		{
			write_version_string(std::cout, config.has("verbose"));
			exit(0);
		}

		if (config.has("help"))
		{
			std::cerr << config << std::endl;
			exit(config.has("help") ? 0 : 1);
		}

		if (config.has("quiet"))
			cif::VERBOSE = -1;
		else
			cif::VERBOSE = config.count("verbose");

		config.parse_config_file("config", "alphafill.conf", { fs::current_path().string(), "/etc/" });

		// --------------------------------------------------------------------

		if (not config.has("pdb-fasta"))
		{
			std::cout << "fasta file not specified" << std::endl;
			exit(1);
		}

		if (not config.has("pdb-dir"))
		{
			std::cout << "PDB directory not specified" << std::endl;
			exit(1);
		}

		// --------------------------------------------------------------------

		std::string fasta = config.get<std::string>("pdb-fasta");
		if (not fs::exists(fasta))
			throw std::runtime_error("PDB-Fasta file does not exist (" + fasta + ")");

		fs::path pdbDir = config.get<std::string>("pdb-dir");
		if (not fs::is_directory(pdbDir))
			throw std::runtime_error("PDB directory does not exist");

		// --------------------------------------------------------------------
		
		if (config.has("validate-fasta"))
			return validateFastA(config.get<std::string>("pdb-fasta"), config.get<std::string>("pdb-dir"), std::thread::hardware_concurrency());

		if (config.has("test-af-id"))
		{
			data_service::instance().exists_in_afdb(config.get<std::string>("test-af-id"));
			return 0;
		}

		// --------------------------------------------------------------------


		std::string command;
		if (not config.operands().empty())
			command = config.operands().front();

		if (command == "server")
			result = server_main(argc - 1, argv + 1);
		else if (command == "process")
			result = alphafill_main(argc - 1, argv + 1);
		else if (command == "validate")
			std::cerr << "unimplemented" << std::endl;
		else if (command == "prepare-pdb-list")
			result = GeneratePDBList();
		else 
		{
			std::cerr << "Usage: alphafill command [options...]" << std::endl
					  << "  where command is one of: 'server', 'process', 'validate', 'prepare-pdb-list'" << std::endl;
			
			exit(1);
		}
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
