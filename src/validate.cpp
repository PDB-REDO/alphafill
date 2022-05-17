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
using Point = mmcif::Point;

// --------------------------------------------------------------------

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

// --------------------------------------------------------------------

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

double Align(const mmcif::Structure &a, mmcif::Structure &b,
	std::vector<Point> &cAlphaA, std::vector<Point> &cAlphaB)
{
	auto ta = CenterPoints(cAlphaA);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << std::endl;

	auto tb = CenterPoints(cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << std::endl;

	auto rotation = AlignPoints(cAlphaB, cAlphaA);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = mmcif::QuaternionToAngleAxis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	b.translateRotateAndTranslate(-tb, rotation, ta);

	for (auto &pt : cAlphaB)
		pt.rotate(rotation);

	double result = CalculateRMSD(cAlphaA, cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << std::endl;

	return result;
}

// --------------------------------------------------------------------
// A blast like alignment. Take two mmcif::Polymers and remove residues
// that are not common in both

std::tuple<std::vector<mmcif::Monomer *>, std::vector<mmcif::Monomer *>> AlignAndTrimSequences(mmcif::Polymer &rx, mmcif::Polymer &ry)
{
	using namespace boost::numeric::ublas;

	int dimX = static_cast<int>(rx.size());
	int dimY = static_cast<int>(ry.size());
	if (dimY == 0 or dimX == 0)
		throw std::runtime_error("Empty chains");

	matrix<float> B(dimX, dimY), Ix(dimX, dimY), Iy(dimX, dimY);
	matrix<int8_t> tb(dimX, dimY);

	int x, y;

	const float
		kMatchReward = 5,
		kMismatchCost = -10,
		kGapOpen = 10, gapExtend = 0.1f;

	float high = 0;
	int highX = 0, highY = 0;

	for (x = 0; x < dimX; ++x)
	{
		for (y = 0; y < dimY; ++y)
		{
			auto &a = rx[x];
			auto &b = ry[y];

			float Ix1 = x > 0 ? Ix(x - 1, y) : 0;
			float Iy1 = y > 0 ? Iy(x, y - 1) : 0;

			// score for alignment
			float M;
			if (a.compoundID() == b.compoundID())
				M = kMatchReward;
			else
				M = kMismatchCost;

			// gap open cost is zero if the PDB ATOM records indicate that a gap
			// should be here.
			float gapOpen = kGapOpen;
			if (y == 0 or (y + 1 < dimY and ry[y + 1].seqID() > ry[y].seqID() + 1))
				gapOpen = 0;

			if (x > 0 and y > 0)
				M += B(x - 1, y - 1);

			float s;
			if (M >= Ix1 and M >= Iy1)
			{
				tb(x, y) = 0;
				B(x, y) = s = M;

				Ix(x, y) = M - (x < dimX - 1 ? gapOpen : 0);
				Iy(x, y) = M - (y < dimY - 1 ? gapOpen : 0);
			}
			else if (Ix1 >= Iy1)
			{
				tb(x, y) = 1;
				B(x, y) = s = Ix1;

				Ix(x, y) = Ix1 - gapExtend;
				Iy(x, y) = M - (y < dimY - 1 ? gapOpen : 0);
				if (Iy(x, y) < Iy1 - gapExtend)
					Iy(x, y) = Iy1 - gapExtend;
			}
			else
			{
				tb(x, y) = -1;
				B(x, y) = s = Iy1;

				Ix(x, y) = M - (x < dimX - 1 ? gapOpen : 0);
				if (Ix(x, y) < Ix1 - gapExtend)
					Ix(x, y) = Ix1 - gapExtend;
				Iy(x, y) = Iy1 - gapExtend;
			}

			if (/*(x == dimX - 1 or y == dimY - 1) and */ high < s)
			{
				high = s;
				highX = x;
				highY = y;
			}
		}
	}

	// assign numbers
	x = highX;
	y = highY;

	std::vector<mmcif::Monomer *> ra, rb;

	while (x >= 0 and y >= 0)
	{
		switch (tb(x, y))
		{
			case -1:
				--y;
				break;

			case 1:
				--x;
				break;

			case 0:
				if (rx[x].compoundID() == ry[y].compoundID())
				{
					ra.emplace_back(&rx[x]);
					rb.emplace_back(&ry[y]);
				}

				--x;
				--y;
		}
	}

	std::reverse(ra.begin(), ra.end());
	std::reverse(rb.begin(), rb.end());

	return {ra, rb};
}

// --------------------------------------------------------------------

std::tuple<std::vector<mmcif::Atom>, std::vector<mmcif::Atom>, std::vector<mmcif::Atom>, std::vector<mmcif::Atom>>
FindAtomsNearLigand(const std::vector<mmcif::Monomer *> &pa, const std::vector<mmcif::Monomer *> &pb,
	const mmcif::Residue &ra, const mmcif::Residue &rb, float maxDistance, const Ligand &ligand)
{
	float maxDistanceSq = maxDistance * maxDistance;

	std::vector<mmcif::Atom> aL, bL, aP, bP;

	std::vector<std::tuple<int, mmcif::Atom>> aI, bI;

	for (auto atom : ra.atoms())
	{
		aL.emplace_back(atom);

		int i = 0;

		for (auto &r : pa)
		{
			for (auto ra : r->atoms())
			{
				if (mmcif::DistanceSquared(atom, ra) <= maxDistanceSq)
					aI.emplace_back(i, ra);
			}

			++i;
		}
	}

	for (auto atom : rb.atoms())
	{
		if (ligand.drops(atom.labelAtomID()))
			continue;

		bL.emplace_back(atom);

		int i = 0;

		for (auto &r : pb)
		{
			for (auto ra : r->atoms())
			{
				if (mmcif::DistanceSquared(atom, ra) <= maxDistanceSq)
					bI.emplace_back(i, ra);
			}

			++i;
		}
	}

	assert(aL.size() == bL.size());
	for (size_t i = 0; i < aL.size(); ++i)
		assert(aL[i].labelAtomID() == ligand.map(bL[i].labelAtomID()));

	auto atomLess = [](const std::tuple<int, mmcif::Atom> &a, const std::tuple<int, mmcif::Atom> &b)
	{
		const auto &[ai, aa] = a;
		const auto &[bi, ba] = b;

		int d = ai - bi;

		if (d == 0)
			d = aa.labelAtomID().compare(ba.labelAtomID());

		return d < 0;
	};

	sort(aI.begin(), aI.end(), atomLess);
	sort(bI.begin(), bI.end(), atomLess);

	auto ai = aI.begin(), bi = bI.begin();

	while (ai != aI.end() or bi != bI.end())
	{
		if (ai == aI.end())
		{
			bi = bI.erase(bi);
			continue;
		}

		if (bi == bI.end())
		{
			ai = aI.erase(ai);
			continue;
		}

		if (atomLess(*ai, *bi))
		{
			ai = aI.erase(ai);
			continue;
		}

		if (atomLess(*bi, *ai))
		{
			bi = bI.erase(bi);
			continue;
		}

		aP.emplace_back(std::get<1>(*ai));
		bP.emplace_back(std::get<1>(*bi));

		++ai;
		++bi;
	}

	return {aP, bP, aL, bL};
}

// --------------------------------------------------------------------

std::tuple<float,float,float> validateCif(cif::Datablock &db, const std::string &asymID, zeep::json::element &info,
	float maxLigandPolyAtomDistance)
{
	std::string pdbID;

	for (auto hit : info["hits"])
	{
		for (auto transplant : hit["transplants"])
		{
			if (transplant["asym_id"].as<std::string>() == asymID)
			{
				pdbID = hit["pdb_id"].as<std::string>();
				break;
			}
		}
	}

	mmcif::File pdbFile(file_locator::get_pdb_file(pdbID));
	mmcif::Structure pdbStructure(pdbFile);
	mmcif::Structure afStructure(db);

	auto &ligands = LigandsTable::instance();

	for (auto hit : info["hits"])
	{
		if (hit["pdb_id"] != pdbID)
			continue;

		for (auto transplant : hit["transplants"])
		{
			std::string asymID = transplant["asym_id"].as<std::string>();

			auto pdbAsymID = hit["pdb_asym_id"].as<std::string>();
			auto pdbCompoundID = transplant["compound_id"].as<std::string>();

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

			return { rmsd1, rmsd2, rmsd3 };
		}
	}

	return { std::nanf("0"), std::nanf("0"), std::nanf("0") };
}
