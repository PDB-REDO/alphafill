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

#include <fstream>
#include <iomanip>

#include <cif++/CifUtils.hpp>
#include <cif++/Structure.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/iostreams/filtering_stream.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
#include <boost/program_options.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "revision.hpp"
#include "utilities.hpp"
// #include "ligands.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

int a_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] uniprot-id [output-file]"s);

	visible_options.add_options()
		("af-id", po::value<std::string>(), "AlphaFold ID");

	po::options_description hidden_options("hidden options");

	po::positional_options_description p;
	p.add("af-id", 1);

	po::variables_map vm = load_options(argc, argv, visible_options, hidden_options, p, "alphafill.conf");

	// --------------------------------------------------------------------

	fs::path dbDir = vm["db-dir"].as<std::string>();
	if (not fs::is_directory(dbDir))
		throw std::runtime_error("AlphfaFill data directory does not exist");

	fs::path pdbDir = vm["pdb-dir"].as<std::string>();
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	file_locator::init(dbDir, vm["structure-name-pattern"].as<std::string>(), vm["metadata-name-pattern"].as<std::string>());

	// --------------------------------------------------------------------

	if (vm.count("af-id") == 0)
	{
		std::cout << "AlphaFold ID not specified" << std::endl;
		exit(1);
	}

	return 0;
}
