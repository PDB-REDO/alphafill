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

#include <boost/numeric/ublas/matrix.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "data-service.hpp"
#include "ligands.hpp"
#include "main.hpp"
#include "revision.hpp"
#include "utilities.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;

using json = zeep::json::element;

// --------------------------------------------------------------------

int validate_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	auto &config = load_and_init_config(
		"alphafill validate [options]",

		mcfp::make_option<std::string>("af-id", "AlphaFold ID"),
		mcfp::make_option<std::string>("pdb-id", "ID of the PDB file"),

		mcfp::make_option<std::string>("db-dir", "Directory containing the alphafilled data"),
		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),

		mcfp::make_option<std::string>("ligands", "af-ligands.cif", "File in CIF format describing the ligands and their modifications"),

		mcfp::make_option<float>("max-ligand-to-polymer-atom-distance,d", 6,
			"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure"),
			
		mcfp::make_hidden_option<std::string>("custom-dir", (fs::temp_directory_path() / "alphafill").string(), "Directory for custom built entries")
			
	);

	parse_argv(argc, argv, config);

	// --------------------------------------------------------------------

	if (not config.has("db-dir"))
	{
		std::cout << "AlphaFill data directory not specified\n";
		return 1;
	}

	if (not config.has("pdb-dir"))
	{
		std::cout << "PDB directory not specified\n";
		return 1;
	}

	// --------------------------------------------------------------------

	if (not config.has("af-id") or not config.has("pdb-id"))
	{
		std::cout << "AlphaFold ID or ligand not specified\n";
		return 1;
	}

	// --------------------------------------------------------------------

	fs::path ligandsFile = config.get("ligands");
	if (not fs::exists(ligandsFile))
	{
		std::cerr << "Ligands file not found\n";
		return 1;
	}

	LigandsTable ligands(ligandsFile);

	// --------------------------------------------------------------------

	float maxLigandPolyAtomDistance = config.get<float>("max-ligand-to-polymer-atom-distance");

	// --------------------------------------------------------------------

	const auto &[type, afID, chunk, version] = parse_af_id(config.get("af-id"));

	cif::file afFile(file_locator::get_structure_file(type, afID, chunk, version));
	cif::mm::structure afStructure(afFile);

	auto pdbID = config.get("pdb-id");
	cif::file pdbFile(file_locator::get_pdb_file(pdbID));
	cif::mm::structure pdbStructure(pdbFile);

	std::ifstream metadata(file_locator::get_metadata_file(afID, chunk, version));

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

				std::vector<cif::point> caA, caP;

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
						std::cerr << "No CA atoms mapped, skipping\n";
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
							std::cerr << ex.what() << '\n';
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
						std::cerr << "Could not find poly atoms near " << afRes << '\n';
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
						  << lA.size() << '\n';
			}
			catch (const std::exception &ex)
			{
				if (cif::VERBOSE >= 0)
					std::cerr << "Failed to process asym " << asymID << " in " << pdbID << '\n'
							  << ex.what() << '\n';

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
						  << 0 << '\n';
			}
		}
	}

	return 0;
}