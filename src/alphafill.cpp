/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

#include <thread>
#include <filesystem>

#include <cif++/Structure.hpp>

#include <boost/program_options.hpp>

#include "blast.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;

// --------------------------------------------------------------------

namespace {
	std::string gVersionNr, gVersionDate;
}

void load_version_info()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)"),
		rxVersionNr2(R"(mkdssp-version: (\d+(?:\.\d+)+))");

#include "revision.hpp"

	struct membuf : public std::streambuf
	{
		membuf(char* data, size_t length)       { this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(kRevision), sizeof(kRevision));

	std::istream is(&buffer);

	std::string line;

	while (getline(is, line))
	{
		std::smatch m;

		if (std::regex_match(line, m, rxVersionNr))
		{
			gVersionNr = m[1];
			if (m[2].matched)
				gVersionNr += '*';
			continue;
		}

		if (std::regex_match(line, m, rxVersionDate))
		{
			gVersionDate = m[1];
			continue;
		}

		// always the first, replace with more specific if followed by the other info
		if (std::regex_match(line, m, rxVersionNr2))
		{
			gVersionNr = m[1];
			continue;
		}
	}
}

std::string get_version_nr()
{
	return gVersionNr/* + '/' + cif::get_version_nr()*/;
}

std::string get_version_date()
{
	return gVersionDate;
}

std::string get_version_string()
{
	return gVersionNr + " " + gVersionDate;
}

// --------------------------------------------------------------------

using mmcif::Point;

std::vector<Point> getCAlphaForChain(const cif::Datablock& db, const std::string& asym_id)
{
	
}

// --------------------------------------------------------------------

int a_main(int argc, const char* argv[])
{
	using namespace std::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);
	visible_options.add_options()
		("dict",		po::value<std::vector<std::string>>(),
															"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")

		("pdb-fasta",	po::value<std::string>(),			"The PDB-REDO fasta file")

		("help,h",											"Display help message")
		("version",											"Print version")
		("verbose,v",										"verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("xyzin,i",				po::value<std::string>(),	"coordinates file")
		("output,o",            po::value<std::string>(),	"Output to this file")
		("debug,d",				po::value<int>(),			"Debug level (for even more verbose output)")

		// ("compounds",			po::value<std::string>(),	"Location of the components.cif file from CCD")
		// ("components",			po::value<std::string>(),	"Location of the components.cif file from CCD, alias")
		// ("extra-compounds",		po::value<std::string>(),	"File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file")
		// ("mmcif-dictionary",	po::value<std::string>(),	"Path to the mmcif_pdbx.dic file to use instead of default")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << get_version_string() << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cerr << visible_options << std::endl;
		exit(0);
	}
	
	if (vm.count("xyzin") == 0)
	{
		std::cerr << "Input file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("pdb-fasta") == 0)
	{
		std::cerr << "PDB-REDO fasta file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("output-format") and vm["output-format"].as<std::string>() != "dssp" and vm["output-format"].as<std::string>() != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// --------------------------------------------------------------------

	// Load extra CCD definitions, if any

	if (vm.count("compounds"))
		cif::addFileResource("components.cif", vm["compounds"].as<std::string>());
	else if (vm.count("components"))
		cif::addFileResource("components.cif", vm["components"].as<std::string>());
	
	if (vm.count("extra-compounds"))
		mmcif::CompoundFactory::instance().pushDictionary(vm["extra-compounds"].as<std::string>());
	
	// And perhaps a private mmcif_pdbx dictionary

	if (vm.count("mmcif-dictionary"))
		cif::addFileResource("mmcif_pdbx_v50.dic", vm["mmcif-dictionary"].as<std::string>());

	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<std::vector<std::string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	// --------------------------------------------------------------------
	
	std::string fasta = vm["pdb-fasta"].as<std::string>();
	if (not fs::exists(fasta))
		throw std::runtime_error("PDB-Fasta file does not exist (" + fasta + ")");

	// --------------------------------------------------------------------
	
	mmcif::File f(vm["xyzin"].as<std::string>());
	// mmcif::Structure structure(f, 1, mmcif::StructureOpenOptions::SkipHydrogen);

	// fetch the (single) chain

	auto& db = f.data();
	for (auto r: db["entity_poly"])
	{
		auto&& [id, seq] = r.get<std::string,std::string>({"entity_id", "pdbx_seq_one_letter_code"});
		
		std::cout << "Blasting:" << std::endl
				  << seq << std::endl
				  << std::endl;

		auto result = BlastP(fasta, seq, "blastp", "BLOSUM62", 0, 10, true, true, -1, -1, 50, std::thread::hardware_concurrency());

		std::cout << "Found " << result.size() << " hits" << std::endl;

		for (auto &hit: result)
			std::cout << hit.id << '\t'
					  << hit.score << '\t'
					  << hit.seq << std::endl;

	}




	return 0;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what (const std::exception& e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception& nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

int main(int argc, const char* argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::addDataDirectory(DATA_DIR);
#endif
		load_version_info();

		result = a_main(argc, argv);
	}
	catch (const std::exception& ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
