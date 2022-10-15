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

#include <cif++.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "revision.hpp"
#include "utilities.hpp"
#include "data-service.hpp"
#include "ligands.hpp"
#include "validate.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
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
	cif::point loc;
	bool ligand;
	std::string compoundID;
	std::string atomID;
};

const cif::mm::residue &guessResidueForLigand(const cif::mm::structure &afs, const std::string &afLigandAsymID, const cif::mm::structure &pdb, const Ligand &ligand)
{
	auto &r1 = afs.get_residue(afLigandAsymID);

	std::string compoundID;
	if (ligand)
		compoundID = ligand.ID();
	if (compoundID.empty())
		compoundID = r1.get_compound_id();

	float sB1 = 0;
	std::vector<cif::point> r1p;

	for (auto a : r1.atoms())
	{
		sB1 += a.get_property_float("B_iso_or_equiv");
		r1p.emplace_back(a.get_location());
	}

	auto c1 = centroid(r1p);

	auto a1 = r1.atoms();

	sort(a1.begin(), a1.end(), [](const cif::mm::atom &a, const cif::mm::atom &b)
		{ return a.get_label_atom_id().compare(b.get_label_atom_id()) < 0; });

	// Key type is difference in sum of b-factors and then distance from centroid
	using M_t = std::tuple<float,float,const cif::mm::residue *>;
	std::vector<M_t> m;
	auto lessM = [](const M_t &a, const M_t &b)
	{
		auto d = std::get<0>(a) - std::get<0>(b);
		if (d == 0)
			d = std::get<1>(a) - std::get<1>(b);
		return d < 0;
	};

	for (auto &r2 : pdb.non_polymers())
	{
		if (r2.get_compound_id() != compoundID)
			continue;

		auto a2 = r2.atoms();

		if (a1.size() != ligand.atom_count(r2))
			continue;

		std::vector<cif::point> r2p;
		float sB2 = 0;

		for (auto a : r2.atoms())
		{
			if (ligand.drops(a.get_label_atom_id()))
				continue;

			sB2 += a.get_property_float("B_iso_or_equiv");
			r2p.emplace_back(a.get_location());
		}

		auto c2 = centroid(r2p);

		m.emplace_back(std::abs(sB2 - sB1), distance(c1, c2), &r2);
		std::push_heap(m.begin(), m.end(), lessM);
	}

	if (m.empty())
		throw std::runtime_error("Could not locate ligand " + afLigandAsymID + " (" + r1.get_compound_id() + ')');

	std::sort_heap(m.begin(), m.end(), lessM);

	return *std::get<2>(m.front());
}

// --------------------------------------------------------------------

using cif::point;
using cif::quaternion;

std::vector<cif::mm::residue *> get_residuesForChain(cif::mm::structure &structure, const std::string &chain_id)
{
	std::vector<cif::mm::residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.get_asym_id() != chain_id)
			continue;

		for (auto &res : poly)
			result.emplace_back(&res);
	}

	return result;
}

std::vector<point> getCAlphaForChain(const std::vector<cif::mm::residue *> &residues)
{
	std::vector<point> result;

	for (auto res : residues)
		result.push_back(res->get_atom_by_atom_id("CA").get_location());

	return result;
}

// --------------------------------------------------------------------

int a_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);

	visible_options.add_options()
		("af-id", po::value<std::string>(), "AlphaFold ID")
		("pdb-id", po::value<std::string>(), "ID of the PDB file")

		("max-ligand-to-polymer-atom-distance,d", po::value<float>()->default_value(6),
			"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure");

	po::options_description hidden_options("hidden options");

	po::positional_options_description p;
	p.add("af-id", 1);
	p.add("pdb-id", 1);

	po::variables_map vm = load_options(argc, argv, visible_options, hidden_options, p, "alphafill.conf");

	// --------------------------------------------------------------------

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

	file_locator::init(vm);

	// --------------------------------------------------------------------

	if (vm.count("af-id") == 0 or vm.count("pdb-id") == 0)
	{
		std::cout << "AlphaFold ID or ligand not specified" << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------

	fs::path ligandsFile = vm["ligands"].as<std::string>();
	if (not fs::exists(ligandsFile))
	{
		std::cerr << "Ligands file not found" << std::endl;
		exit(1);
	}

	LigandsTable ligands(ligandsFile);

	// --------------------------------------------------------------------
	
	float maxLigandPolyAtomDistance = vm["max-ligand-to-polymer-atom-distance"].as<float>();

	// --------------------------------------------------------------------

	const auto &[afID, chunk] = parse_af_id(vm["af-id"].as<std::string>());

	cif::file afFile(file_locator::get_structure_file(afID, chunk));
	cif::mm::structure afStructure(afFile);

	auto pdbID = vm["pdb-id"].as<std::string>();
	cif::file pdbFile(file_locator::get_pdb_file(pdbID));
	cif::mm::structure pdbStructure(pdbFile);

	std::ifstream metadata(file_locator::get_metdata_file(afID, chunk));

	json info;
	zeep::json::parse_json(metadata, info);

	for (auto hit : info["hits"])
	{
		if (hit["pdb_id"] != pdbID)
			continue;

		for (auto transplant : hit["transplants"])
		{
			std::string asymID = transplant["asym_id"].as<std::string>();

			auto pdbAsymID = hit["pdb_asym_id"].as<std::string>();
			auto pdbCompoundID = transplant["compound_id"].as<std::string>();

			try
			{
				auto ligand = ligands[pdbCompoundID];

				auto &afPolyS = afStructure.get_polymer_by_asym_id("A");
				auto &pdbPolyS = pdbStructure.get_polymer_by_asym_id(pdbAsymID);

				auto &&[afPoly, pdbPoly] = AlignAndTrimSequences(afPolyS, pdbPolyS);

				if (afPoly.size() != pdbPoly.size())
					throw std::runtime_error("polymers differ in length");

				std::vector<point> caA, caP;

				for (size_t i = 0; i < afPoly.size(); ++i)
				{
					auto af_ca = afPoly[i]->get_atom_by_atom_id("CA");
					if (not af_ca)
						continue;

					auto pdb_ca = pdbPoly[i]->get_atom_by_atom_id("CA");
					if (not pdb_ca)
						continue;

					caA.push_back(af_ca.get_location());
					caP.push_back(pdb_ca.get_location());
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
					try
					{
						auto &rA = *afPoly[i];
						auto &rP = *pdbPoly[i];

						if (rA.get_compound_id() == "TYR" or rA.get_compound_id() == "PHE")
						{
							if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
							{
								pdbStructure.swap_atoms(rP.get_atom_by_atom_id("CD1"), rP.get_atom_by_atom_id("CD2"));
								pdbStructure.swap_atoms(rP.get_atom_by_atom_id("CE1"), rP.get_atom_by_atom_id("CE2"));
							}

							continue;
						}

						if (rA.get_compound_id() == "ASP")
						{
							if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
								pdbStructure.swap_atoms(rP.get_atom_by_atom_id("OD1"), rP.get_atom_by_atom_id("OD2"));
							continue;
						}

						if (rA.get_compound_id() == "GLU")
						{
							if (std::abs(rA.chi(2) - rP.chi(2)) > 90)
								pdbStructure.swap_atoms(rP.get_atom_by_atom_id("OE1"), rP.get_atom_by_atom_id("OE2"));
							continue;
						}
					}
					catch (const std::exception &ex)
					{
						if (cif::VERBOSE > 0)
							std::cerr << ex.what() << std::endl;
					}
				}

				// Align the PDB structure on the AF structure, based on C-alpha
				Align(afStructure, pdbStructure, caA, caP);

				auto &afRes = afStructure.get_residue(asymID);
				// auto &pdbRes = guessResidueForLigand(afStructure, asymID, pdbStructure, ligand);
				auto &pdbRes = pdbStructure.get_residue(transplant["pdb_asym_id"].as<std::string>());

				if (afRes.get_compound_id() != pdbRes.get_compound_id())
					throw std::runtime_error("Compound ID's do not match: " + afRes.get_compound_id() + " != " + pdbRes.get_compound_id());

				// collect atoms around ligand

				auto &&[pA, pP, lA, lP] = FindAtomsNearLigand(afPoly, pdbPoly, afRes, pdbRes, maxLigandPolyAtomDistance, ligand);

				if (pA.empty())
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Could not find poly atoms near " << afRes << std::endl;
					continue;
				}

				std::vector<cif::mm::atom> cA = pA, cP = pP;
				cA.insert(cA.end(), lA.begin(), lA.end());
				cP.insert(cP.end(), lP.begin(), lP.end());

				auto rmsd1 = Align(cA, cP);

				auto rmsd2 = CalculateRMSD(pA, pP);
				auto rmsd3 = CalculateRMSD(lA, lP);

				std::cout << afID << '\t'
						<< pdbID << '\t'
						<< afRes.get_compound_id() << '\t'
						<< pdbRes.get_auth_asym_id() << '\t'
						<< pdbRes.get_auth_seq_id() << '\t'
						<< pdbRes.get_pdb_ins_code() << '\t'
						<< std::setprecision(5) << rmsd1 << '\t'
						<< std::setprecision(5) << rmsd2 << '\t'
						<< std::setprecision(5) << rmsd3 << '\t'
						<< pA.size() << '\t'
						<< lA.size() << std::endl;
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE >= 0)
					std::cerr << "Failed to process asym " << asymID << " in " << pdbID << std::endl
							  << ex.what() << std::endl;

				std::cout << afID << '\t'
						<< pdbID << '\t'
						<< '"' << transplant["compound_id"].as<std::string>() << '/' << transplant["analogue_id"].as<std::string>() << '"' << '\t'
						<< "\"?\"" << '\t'
						<< 0 << '\t'
						<< "\"?\"" << '\t'
						<< 0 << '\t'
						<< 0 << '\t'
						<< 0 << '\t'
						<< 0 << '\t'
						<< 0 << std::endl;
			}
		}
	}

	return 0;
}