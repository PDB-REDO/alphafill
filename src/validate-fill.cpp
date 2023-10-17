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

std::tuple<double,double,double,size_t,size_t>
ValidateTransplant(const LigandsTable &ligands, std::string_view af_id, std::string_view pdb_id,
	std::string_view transpant_auth_asym_id, int transplant_auth_seq_id,
	float maxLigandPolyAtomDistance)
{
	const auto &[type, afID, chunk, version] = parse_af_id(std::string{ af_id });

	std::ifstream metadata(file_locator::get_metadata_file(afID, chunk, version));
	json info;
	zeep::json::parse_json(metadata, info);

	std::string asymID, pdbAsymID, pdbCompoundID;

	for (zeep::json::element &hit : info["hits"])
	{
		if (hit["pdb_id"].as<std::string>() != pdb_id)
			continue;
		
		for (auto transplant : hit["transplants"])
		{
			if (transplant["pdb_auth_asym_id"].as<std::string>() != transpant_auth_asym_id or
				transplant["pdb_auth_seq_id"].as<int>() != transplant_auth_seq_id)
			{
				continue;
			}

			asymID = transplant["asym_id"].as<std::string>();
			pdbAsymID = transplant["pdb_asym_id"].as<std::string>();
			pdbCompoundID = transplant["compound_id"].as<std::string>();

			break;
		}
	}

	if (asymID.empty())
	{
		std::ostringstream os;
		os << "Transplant not found for " << pdb_id << '/'
		   << transpant_auth_asym_id << '/'
		   << transplant_auth_seq_id;

		throw std::runtime_error(os.str());
	}

	cif::file afFile(file_locator::get_structure_file(type, afID, chunk, version));
	cif::mm::structure afStructure(afFile);

	cif::file pdbFile(file_locator::get_pdb_file(std::string{ pdb_id }));
	cif::mm::structure pdbStructure(pdbFile);

	auto &afRes = afStructure.get_residue(asymID);
	auto &pdbRes = pdbStructure.get_residue(pdbAsymID);

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

		return {};
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

	// if (afRes.get_compound_id() != pdbRes.get_compound_id())
	// 	throw std::runtime_error("Compound ID's do not match: " + afRes.get_compound_id() + " != " + pdbRes.get_compound_id());

	// collect atoms around ligand

	auto &&[pA, pP, lA, lP] = FindAtomsNearLigand(afPoly, pdbPoly, afRes, pdbRes, maxLigandPolyAtomDistance, ligand);

	if (pA.empty())
	{
		if (cif::VERBOSE > 0)
			std::cerr << "Could not find poly atoms near " << afRes << '\n';
		return {};
	}

	std::vector<cif::mm::atom> cA = pA, cP = pP;
	cA.insert(cA.end(), lA.begin(), lA.end());
	cP.insert(cP.end(), lP.begin(), lP.end());

	auto rmsd1 = Align(cA, cP);

	auto rmsd2 = CalculateRMSD(pA, pP);
	auto rmsd3 = CalculateRMSD(lA, lP);

	return { rmsd1, rmsd2, rmsd3, pA.size(), lA.size() };
}

int validate_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	auto &config = load_and_init_config(
		"alphafill validate [options] <inputfile> [<outputfile>]",

		// mcfp::make_option<std::string>("af-id", "AlphaFold ID"),
		// mcfp::make_option<std::string>("pdb-id", "ID of the PDB file"),

		mcfp::make_option<std::string>("db-dir", "Directory containing the alphafilled data"),
		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),

		mcfp::make_option<std::string>("ligands", "af-ligands.cif", "File in CIF format describing the ligands and their modifications"),

		mcfp::make_option<float>("max-ligand-to-polymer-atom-distance,d", 6,
			"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure"),

		mcfp::make_option<std::string>("structure-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.gz", "Pattern for locating structure files"),
		mcfp::make_option<std::string>("metadata-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.json", "Pattern for locating metadata files"),
		mcfp::make_option<std::string>("pdb-name-pattern", "${pdb-dir}/${id:1:2}/${id}/${id}_final.cif", "Pattern for locating PDB files"),

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

	// if (not config.has("af-id") or not config.has("pdb-id"))
	if (config.operands().empty())
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

	std::ifstream in(config.operands().front());
	std::ofstream outf;

	std::streambuf *out;
	if (config.operands().size() == 1)
		out = std::cout.rdbuf();
	else
	{
		outf.open(config.operands()[1]);
		if (not outf.is_open())
			throw std::runtime_error("Could not open output file");
		
		out = outf.rdbuf();
	}

	std::ostream s_out(out);

	if (not in.is_open())
		throw std::runtime_error("Could not open input file");

	std::string header;
	getline(in, header);

	auto headers = cif::split(header, ",");
	auto h_ix_1 = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), "AFill_ID"));
	auto h_ix_2 = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), "pdb_id"));
	auto h_ix_3 = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), "pdb_auth_asym_id"));
	auto h_ix_4 = std::distance(headers.begin(), std::find(headers.begin(), headers.end(), "pdb_auth_seq_id"));

	headers.emplace_back("rmsd");
	headers.emplace_back("rmsd-poly-atoms");
	headers.emplace_back("rmsd-ligand-atoms");
	headers.emplace_back("rmsd-poly-atom-count");
	headers.emplace_back("rmsd-ligand-atom-count");

	std::string line;
	while (getline(in, line))
	{
		try
		{
			auto flds = cif::split(line, ",");
			if (flds.empty())
				continue;
			
			auto af_id = flds.at(h_ix_1);
			auto pdb_id = flds.at(h_ix_2);
			auto transplant_auth_asym_id = flds.at(h_ix_3);
			int transplant_auth_seq_id = std::stoi(std::string{ flds.at(h_ix_4) });

			auto r = ValidateTransplant(ligands, af_id, pdb_id, transplant_auth_asym_id, transplant_auth_seq_id, maxLigandPolyAtomDistance);

			s_out << line << '\t'
				  << std::get<0>(r) << '\t'
				  << std::get<1>(r) << '\t'
				  << std::get<2>(r) << '\t'
				  << std::get<3>(r) << '\t'
				  << std::get<4>(r) << '\n';
		}
		catch (const std::exception &ex)
		{
			std::cerr << "Error processing line:\n  " << line << '\n' << ex.what() << std::endl;
		}
	}

	return 0;
}