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
#include <zeep/json/parser.hpp>

#include "revision.hpp"
#include "utilities.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

fs::path pdbFileForID(const fs::path &pdbDir, std::string pdb_id)
{
	for (auto &ch : pdb_id)
		ch = std::tolower(ch);

	// try a PDB-REDO layout first
	fs::path pdb_path = pdbDir / pdb_id.substr(1, 2) / pdb_id / (pdb_id + "_final.cif");
	if (not fs::exists(pdb_path))
		pdb_path = pdbDir / pdb_id.substr(1, 2) / (pdb_id + ".cif.gz");

	if (not fs::exists(pdb_path))
		throw std::runtime_error("PDB file for " + pdb_id + " not found");

	return pdb_path;
}

// --------------------------------------------------------------------

struct VAtom
{
	mmcif::Point loc;
	bool ligand;
	std::string compoundID;
	std::string atomID;
};

const mmcif::Residue &guessResidueForLigand(const mmcif::Structure &afs, const std::string &afLigandAsymID, const mmcif::Structure &pdb)
{
	auto &r1 = afs.getResidue(afLigandAsymID);

	auto a1 = r1.atoms();
	sort(a1.begin(), a1.end(), [](const mmcif::Atom &a, const mmcif::Atom &b)
	{
		return a.labelAtomID().compare(b.labelAtomID());
	});

	for (auto &r2 : pdb.getNonPolymers())
	{
		if (r1.compoundID() != r2.compoundID())
			continue;

		auto a2 = r2.atoms();

		if (a1.size() != a2.size())
			continue;

		sort(a2.begin(), a2.end(), [](const mmcif::Atom &a, const mmcif::Atom &b)
		{
			return a.labelAtomID().compare(b.labelAtomID());
		});

		bool isSame = true;

		for (size_t i = 0; isSame and i < a1.size(); ++i)
		{
			auto &aa1 = a1[i];
			auto &aa2 = a2[i];

			if (aa1.labelAtomID() != aa2.labelAtomID() or
				std::fabs(aa1.get_property<float>("B_iso_or_equiv") - aa2.get_property<float>("B_iso_or_equiv")) >= 0.01f)
			{
				isSame = false;
				break;
			}
		}

		if (isSame)
			return r2;
	}

	throw std::runtime_error("Could not locate ligand " + afLigandAsymID);
}


// --------------------------------------------------------------------

using mmcif::Point;
using mmcif::Quaternion;

std::vector<mmcif::Residue *> getResiduesForChain(mmcif::Structure &structure, const std::string &chain_id)
{
	std::vector<mmcif::Residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.chainID() != chain_id)
			continue;

		for (auto &res : poly)
			result.emplace_back(&res);
	}

	return result;
}

std::vector<Point> getCAlphaForChain(const std::vector<mmcif::Residue *> &residues)
{
	std::vector<Point> result;

	for (auto res : residues)
		result.push_back(res->atomByID("CA").location());

	return result;
}

double CalculateRMSD(const std::vector<mmcif::Point> &pa, const std::vector<mmcif::Point> &pb)
{
	return RMSd(pa, pb);
}

double CalculateRMSD(const std::vector<mmcif::Atom> &a, const std::vector<mmcif::Atom> &b)
{
	std::vector<Point> pa, pb;

	for (auto &atom : a)
		pa.emplace_back(atom.location());

	for (auto &atom : b)
		pb.emplace_back(atom.location());

	return RMSd(pa, pb);
}

double Align(std::vector<mmcif::Atom> &aA, std::vector<mmcif::Atom> &aB)
{
	std::vector<Point> pA, pB;

	for (auto &a : aA)
		pA.emplace_back(a.location());

	for (auto &b : aB)
		pB.emplace_back(b.location());

	auto ta = CenterPoints(pA);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << std::endl;

	auto tb = CenterPoints(pB);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << std::endl;

	auto rotation = AlignPoints(pB, pA);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = mmcif::QuaternionToAngleAxis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	for (auto b : aB)
		b.translateRotateAndTranslate(-tb, rotation, ta);

	for (auto &pt : pB)
		pt.rotate(rotation);

	double result = CalculateRMSD(pA, pB);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << std::endl;

	return result;
}

std::tuple<std::vector<Point>, std::vector<Point>> selectAtomsNearResidue(
	const std::vector<mmcif::Residue *> &pdb, const std::vector<size_t> &pdb_ix,
	const std::vector<mmcif::Residue *> &af, const std::vector<size_t> &af_ix,
	const std::vector<mmcif::Atom> &residue, float maxDistance)
{
	std::vector<Point> ra, rb;

	assert(pdb_ix.size() == af_ix.size());

	for (size_t i = 0; i < pdb_ix.size(); ++i)
	{
		bool nearby = false;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(pdb_ix[i] < pdb.size());

			auto atom = pdb[pdb_ix[i]]->atomByID(atom_id);
			if (not atom)
				continue;

			for (auto &b : residue)
			{
				if (Distance(atom, b) <= maxDistance)
				{
					nearby = true;
					break;
				}
			}

			if (nearby)
				break;
		}

		if (not nearby)
			continue;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(af_ix[i] < af.size());
			assert(pdb_ix[i] < pdb.size());

			auto pt_a = pdb[pdb_ix[i]]->atomByID(atom_id);
			auto pt_b = af[af_ix[i]]->atomByID(atom_id);

			if (not pt_a and pt_b)
				continue;

			ra.push_back(pt_a.location());
			rb.push_back(pt_b.location());
		}
	}

	return {ra, rb};
}

// --------------------------------------------------------------------

std::tuple<std::vector<mmcif::Atom>,std::vector<mmcif::Atom>,std::vector<mmcif::Atom>,std::vector<mmcif::Atom>>
FindAtomsNearLigand(const mmcif::Polymer &pa, const mmcif::Polymer &pb,
	const mmcif::Residue &ra, const mmcif::Residue &rb, float maxDistance)
{
	float dsq = maxDistance * maxDistance;

	std::vector<mmcif::Atom> aL, bL, aP, bP;

	for (auto atom : ra.atoms())
	{
		aL.emplace_back(atom);

		for (auto &r : pa)
		{
			for (auto ra : r.atoms())
			{
				if (mmcif::DistanceSquared(atom, ra) <= dsq)
					aP.emplace_back(ra);
			}
		}
	}
	
	for (auto atom : rb.atoms())
	{
		bL.emplace_back(atom);

		for (auto &r : pb)
		{
			for (auto ra : r.atoms())
			{
				if (mmcif::DistanceSquared(atom, ra) <= dsq)
					bP.emplace_back(ra);
			}
		}
	}

	auto atomLess = [](const mmcif::Atom &a, const mmcif::Atom &b)
	{
		int d = a.labelSeqID() - b.labelSeqID();

		if (d == 0)
			d = a.labelAtomID().compare(b.labelAtomID());
		
		return d < 0;
	};

	sort(aP.begin(), aP.end(), atomLess);
	sort(bP.begin(), bP.end(), atomLess);

	auto ai = aP.begin(), bi = bP.begin();

	while (ai != aP.end() or bi != bP.end())
	{
		if (ai == aP.end())
		{
			bi = bP.erase(bi);
			continue;
		}

		if (bi == bP.end())
		{
			ai = aP.erase(ai);
			continue;
		}

		if (atomLess(*ai, *bi))
		{
			ai = aP.erase(ai);
			continue;
		}

		if (atomLess(*bi, *ai))
		{
			bi = bP.erase(bi);
			continue;
		}

		++ai;
		++bi;
	}

	return { aP, bP, aL, bL };
}

// --------------------------------------------------------------------

int a_main(int argc, const char *argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);

	visible_options.add_options()
		("db-dir",	po::value<std::string>(),	"Directory containing the af-filled data")
		("pdb-dir",	po::value<std::string>(),	"Directory containing the mmCIF files for the PDB")

		("structure-name-pattern",	po::value<std::string>(),	"Pattern for locating structure files")
		("metadata-name-pattern",	po::value<std::string>(),	"Pattern for locating metadata files")

		("af-id",	po::value<std::string>(),	"AlphaFold ID")
		("pdb-id", po::value<std::string>(),	"ID of the PDB file")

		("max-ligand-to-polymer-atom-distance,d", po::value<float>()->default_value(6),
												"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure")

		// ("threads,t", po::value<size_t>()->default_value(1), "Number of threads to use, zero means all available cores")

		("config", po::value<std::string>(), "Config file")
		("help,h", "Display help message")
		("version", "Print version")
		("verbose,v", "Verbose output")
		("quiet", "Do not produce warnings");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d", po::value<int>(), "Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("af-id", 1);
	p.add("pdb-id", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv)
		.options(cmdline_options)
		.positional(p)
		.run(), vm);

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

	fs::path dbDir = vm["db-dir"].as<std::string>();
	if (not fs::is_directory(dbDir))
		throw std::runtime_error("AlphfaFill data directory does not exist");

	fs::path pdbDir = vm["pdb-dir"].as<std::string>();
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	file_locator::init(dbDir, vm["structure-name-pattern"].as<std::string>(), vm["metadata-name-pattern"].as<std::string>());

	// --------------------------------------------------------------------
	
	if (vm.count("af-id") == 0 or vm.count("pdb-id") == 0)
	{
		std::cout << "AlphaFold ID or ligand not specified" << std::endl;
		exit(1);
	}

	float maxLigandPolyAtomDistance = vm["max-ligand-to-polymer-atom-distance"].as<float>();

	// --------------------------------------------------------------------

	auto afID = vm["af-id"].as<std::string>();
	mmcif::File afFile(file_locator::get_structure_file(afID));
	mmcif::Structure afStructure(afFile);

	auto pdbID = vm["pdb-id"].as<std::string>();
	mmcif::File pdbFile(pdbFileForID(pdbDir, pdbID));
	mmcif::Structure pdbStructure(pdbFile);

	std::ifstream metadata(file_locator::get_metdata_file(afID));

	json info;
	zeep::json::parse_json(metadata, info);

	for (auto hit : info["hits"])
	{
		if (hit["pdb_id"] != pdbID)
			continue;

		for (auto transplant : hit["transplants"])
		{
			std::string asymID = transplant["asym_id"].as<std::string>();

			auto &afRes = afStructure.getResidue(asymID);
			auto &pdbRes = guessResidueForLigand(afStructure, asymID, pdbStructure);

			assert(afRes.compoundID() == pdbRes.compoundID());

			auto pdbAsymID = hit["pdb_asym_id"].as<std::string>();

			auto &afPoly = afStructure.getPolymerByAsymID("A");
			auto &pdbPoly = pdbStructure.getPolymerByAsymID(pdbAsymID);

			if (afPoly.size() != pdbPoly.size())
				throw std::runtime_error("polymers differ in length");
			
			std::vector<Point> caA, caP;

			for (size_t i = 0; i < afPoly.size(); ++i)
			{
				auto af_ca = afPoly[i].atomByID("CA");
				if (not af_ca)
					continue;

				auto pdb_ca = pdbPoly[i].atomByID("CA");
				if (not pdb_ca)
					continue;

				caA.push_back(af_ca.location());
				caP.push_back(pdb_ca.location());
			}

			if (caA.empty())
			{
				if (cif::VERBOSE > 0)
					std::cerr << "No CA atoms mapped, skipping" << std::endl;
				continue;
			}			

			// [11-4 11:43] Robbie Joosten
			// Bij de alignment moeten we rekening houden met de conformatie van TYR, PHE, ASP, en GLU.
			// Het verschil in de laatste torsiehoek moet geminimaliseerd worden door de zijketens (in
			// het PDB_REDO model) 180 graden te flippen (i.e. de atoomnamen te swappen). Voor TYR, PHE
			// en ASP gaat het om de torsiehoek chi-2. In GLU gaat het om chi-3.

			for (size_t i = 0; i < afPoly.size(); ++i)
			{
				auto &rA = afPoly[i];
				auto &rP = pdbPoly[i];

				if (rA.compoundID() == "TYR" or rA.compoundID() == "PHE")
				{
					if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
					{
						pdbStructure.swapAtoms(rP.atomByID("CD1"), rP.atomByID("CD2"));
						pdbStructure.swapAtoms(rP.atomByID("CE1"), rP.atomByID("CE2"));
					}

					continue;
				}

				if (rA.compoundID() == "ASP")
				{
					if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
						pdbStructure.swapAtoms(rP.atomByID("OD1"), rP.atomByID("OD2"));
					continue;
				}

				if (rA.compoundID() == "GLU")
				{
					if (std::abs(rA.chi(2) - rP.chi(2)) > 90)
						pdbStructure.swapAtoms(rP.atomByID("OE1"), rP.atomByID("OE2"));
					continue;
				}
			}

			// collect atoms around ligand

			auto &&[pA, pP, lA, lP] = FindAtomsNearLigand(afPoly, pdbPoly, afRes, pdbRes, maxLigandPolyAtomDistance);

			std::vector<mmcif::Atom> cA = pA, cP = pP;
			cA.insert(cA.end(), lA.begin(), lA.end());
			cP.insert(cP.end(), lP.begin(), lP.end());

			auto rmsd1 = Align(cA, cP);

			auto rmsd2 = CalculateRMSD(pA, pP);
			auto rmsd3 = CalculateRMSD(lA, lP);

			std::cout << pdbID << '\t'
					  << pdbAsymID << '\t'
					  << afRes.compoundID() << '\t'
					  << afRes.asymID() << '\t'
					  << afRes.seqID() << '\t'
					  << afRes.authAsymID() << '\t'
					  << afRes.authSeqID() << '\t'
					  << afRes.authInsCode() << '\t'
					  << rmsd1 << '\t'
					  << rmsd2 << '\t'
					  << rmsd3 << '\t'
					  << std::endl;
		}
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
