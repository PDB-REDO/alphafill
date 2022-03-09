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

#include <fstream>
#include <iomanip>

#include <cif++/CifUtils.hpp>
#include <cif++/Structure.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>

#include <zeep/json/element.hpp>

#include "blast.hpp"
#include "queue.hpp"
#include "revision.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

using mmcif::Point;

struct CAtom
{
	CAtom(const CAtom &) = default;
	CAtom(CAtom &&) = default;

	CAtom(mmcif::AtomType type, Point pt, int charge)
		: type(type), pt(pt)
	{
		radius = charge == 0 ?
			mmcif::AtomTypeTraits(type).radius(mmcif::RadiusType::VanderWaals) : 
			mmcif::AtomTypeTraits(type).effective_ionic_radius(charge);

		if (std::isnan(radius))
			throw std::runtime_error("Unknown radius for atom " + mmcif::AtomTypeTraits(type).symbol() + " with charge " + std::to_string(charge));
	}

	CAtom(const mmcif::Atom &atom)
		: CAtom(atom.type(), atom.location(), atom.charge())
	{
	}

	mmcif::AtomType type;
	Point pt;
	float radius;
};

// --------------------------------------------------------------------

json CalculateClashScore(const std::vector<CAtom> &polyAtoms, const std::vector<CAtom> &resAtoms, float maxDistance)
{
	auto maxDistanceSq = maxDistance * maxDistance;

	json result;

	auto &dp = result["distance"];

	int n = 0, o = 0;

	for (auto &pa : polyAtoms)
	{
		bool near = false;

		for (auto &ra : resAtoms)
		{
			auto d = DistanceSquared(pa.pt, ra.pt);

			if (d >= 2 * maxDistanceSq)
				continue;
			
			near = true;

			d = std::sqrt(d);

			auto overlap = pa.radius + ra.radius - d;
			if (overlap < 0)
				overlap = 0;
			
			if (overlap > 0)
				++n;
			
			json d_pair;
			d_pair.push_back(d);
			d_pair.push_back(overlap);

			dp.push_back(std::move(d_pair));
		}

		if (near)
			++o;
	}

	result["clash_count"] = n;
	result["poly_atom_count"] = o;
	result["ligand_atom_count"] = resAtoms.size();

	return result;
}

// --------------------------------------------------------------------

int a_main(int argc, const char *argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);

	visible_options.add_options()
		("distance-cutoff", po::value<float>()->default_value(4), "The max distance between polymer atoms and ligand atoms used in calculating clash scores")
		("config", po::value<std::string>(), "Config file")
		("help,h", "Display help message")
		("version", "Print version")
		("verbose,v", "Verbose output")
		("quiet", "Do not produce warnings");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("xyzin,i", po::value<std::string>(), "coordinates file")
		("output,o", po::value<std::string>(), "Output to this file")
		("debug,d", po::value<int>(), "Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "alphafill.conf";
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

	float maxDistance = vm["distance-cutoff"].as<float>();

	// --------------------------------------------------------------------
	
	fs::path xyzin = vm["xyzin"].as<std::string>();
	mmcif::File f(xyzin);
	mmcif::Structure structure(f, 1, mmcif::StructureOpenOptions::SkipHydrogen);
	
	auto &db = f.data();

	auto &polies = structure.polymers();
	if (polies.size() != 1)
		throw std::runtime_error("Number of polymers in this file is not exactly 1");

	auto &poly = polies.front();

	std::vector<CAtom> polyAtoms;

	for (auto &res : poly)
	{
		for (auto &atom : res.atoms())
			polyAtoms.emplace_back(atom);
	}

	std::vector<std::string> ligandIDs;
	for (auto &&[id] : db["struct_asym"].rows<std::string>("id"))
	{
		if (id == poly.asymID())
			continue;
		ligandIDs.emplace_back(id);
	}

	for (auto &asymID : ligandIDs)
	{
		std::vector<CAtom> resAtoms;

		// auto &res = structure.getResidue(asymID);

		// for (auto &atom : res.atoms())
		// 	resAtoms.emplace_back(atom);

		// assert(not resAtoms.empty());

		for (const auto &[type_symbol, x, y, z, charge, comp_id] :
			db["atom_site"].find<std::string,float,float,float,std::optional<int>,std::string>(
				"label_asym_id"_key == asymID,
				"type_symbol", "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_formal_charge", "label_comp_id"))
		{
			mmcif::AtomTypeTraits att(type_symbol);

			int formal_charge = charge.value_or(0);

			if (not charge.has_value() and att.isMetal())
			{
				auto compound = mmcif::CompoundFactory::instance().create(comp_id);
				if (compound)
					formal_charge = compound->formalCharge();
			}

			resAtoms.emplace_back(att.type(), Point{ x, y, z }, formal_charge);
		}

		json clash{
			{ "asym_id", asymID },
			{ "clash", CalculateClashScore(polyAtoms, resAtoms, maxDistance) }
		};

		std::cout << clash << std::endl;
	}

	return 0;
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

int main(int argc, const char *argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::addDataDirectory(DATA_DIR);
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
