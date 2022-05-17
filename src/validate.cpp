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
#include <boost/numeric/ublas/matrix.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "utilities.hpp"
#include "data-service.hpp"
#include "ligands.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

std::tuple<float,float,float> validateCif(const cif::Datablock &db, zeep::json::element &info)
{
	mmcif::File pdbFile(file_locator::get_pdb_file(pdbID));
	mmcif::Structure pdbStructure(pdbFile);


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

				auto &afPolyS = afStructure.getPolymerByAsymID("A");
				auto &pdbPolyS = pdbStructure.getPolymerByAsymID(pdbAsymID);

				auto &&[afPoly, pdbPoly] = AlignAndTrimSequences(afPolyS, pdbPolyS);

				if (afPoly.size() != pdbPoly.size())
					throw std::runtime_error("polymers differ in length");

				std::vector<Point> caA, caP;

				for (size_t i = 0; i < afPoly.size(); ++i)
				{
					auto af_ca = afPoly[i]->atomByID("CA");
					if (not af_ca)
						continue;

					auto pdb_ca = pdbPoly[i]->atomByID("CA");
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
					try
					{
						auto &rA = *afPoly[i];
						auto &rP = *pdbPoly[i];

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
					catch (const std::exception &ex)
					{
						if (cif::VERBOSE > 0)
							std::cerr << ex.what() << std::endl;
					}
				}

				// Align the PDB structure on the AF structure, based on C-alpha
				Align(afStructure, pdbStructure, caA, caP);

				auto &afRes = afStructure.getResidue(asymID);
				// auto &pdbRes = guessResidueForLigand(afStructure, asymID, pdbStructure, ligand);
				auto &pdbRes = pdbStructure.getResidue(transplant["pdb_asym_id"].as<std::string>());

				if (afRes.compoundID() != pdbRes.compoundID())
					throw std::runtime_error("Compound ID's do not match: " + afRes.compoundID() + " != " + pdbRes.compoundID());

				// collect atoms around ligand

				auto &&[pA, pP, lA, lP] = FindAtomsNearLigand(afPoly, pdbPoly, afRes, pdbRes, maxLigandPolyAtomDistance, ligand);

				if (pA.empty())
				{
					if (cif::VERBOSE > 0)
						std::cerr << "Could not find poly atoms near " << afRes << std::endl;
					continue;
				}

				std::vector<mmcif::Atom> cA = pA, cP = pP;
				cA.insert(cA.end(), lA.begin(), lA.end());
				cP.insert(cP.end(), lP.begin(), lP.end());

				auto rmsd1 = Align(cA, cP);

				auto rmsd2 = CalculateRMSD(pA, pP);
				auto rmsd3 = CalculateRMSD(lA, lP);

				std::cout << afID << '\t'
						<< pdbID << '\t'
						<< afRes.compoundID() << '\t'
						<< pdbRes.authAsymID() << '\t'
						<< pdbRes.authSeqID() << '\t'
						<< pdbRes.authInsCode() << '\t'
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
