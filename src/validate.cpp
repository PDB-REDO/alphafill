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

#include "utilities.hpp"
#include "data-service.hpp"
#include "ligands.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;
namespace ba = boost::algorithm;

using json = zeep::json::element;
using point = cif::point;

// --------------------------------------------------------------------

double CalculateRMSD(const std::vector<cif::point> &pa, const std::vector<cif::point> &pb)
{
	return RMSd(pa, pb);
}

double CalculateRMSD(const std::vector<cif::mm::atom> &a, const std::vector<cif::mm::atom> &b)
{
	std::vector<point> pa, pb;

	for (auto &atom : a)
		pa.emplace_back(atom.get_location());

	for (auto &atom : b)
		pb.emplace_back(atom.get_location());

	return RMSd(pa, pb);
}

// --------------------------------------------------------------------

double Align(std::vector<cif::mm::atom> &aA, std::vector<cif::mm::atom> &aB)
{
	std::vector<point> pA, pB;

	for (auto &a : aA)
		pA.emplace_back(a.get_location());

	for (auto &b : aB)
		pB.emplace_back(b.get_location());

	auto ta = center_points(pA);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << std::endl;

	auto tb = center_points(pB);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << std::endl;

	auto rotation = align_points(pB, pA);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = cif::quaternion_to_angle_axis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	for (auto b : aB)
		b.translate_rotate_and_translate(-tb, rotation, ta);

	for (auto &pt : pB)
		pt.rotate(rotation);

	double result = CalculateRMSD(pA, pB);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << std::endl;

	return result;
}

double Align(const cif::mm::structure &a, cif::mm::structure &b,
	std::vector<point> &cAlphaA, std::vector<point> &cAlphaB)
{
	auto ta = center_points(cAlphaA);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << std::endl;

	auto tb = center_points(cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << std::endl;

	auto rotation = align_points(cAlphaB, cAlphaA);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = cif::quaternion_to_angle_axis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	b.translate_rotate_and_translate(-tb, rotation, ta);

	for (auto &pt : cAlphaB)
		pt.rotate(rotation);

	double result = CalculateRMSD(cAlphaA, cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << std::endl;

	return result;
}

// --------------------------------------------------------------------
// A blast like alignment. Take two cif::mm::polymers and remove residues
// that are not common in both

std::tuple<std::vector<cif::mm::monomer *>, std::vector<cif::mm::monomer *>> AlignAndTrimSequences(cif::mm::polymer &rx, cif::mm::polymer &ry)
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
			if (a.get_compound_id() == b.get_compound_id())
				M = kMatchReward;
			else
				M = kMismatchCost;

			// gap open cost is zero if the PDB ATOM records indicate that a gap
			// should be here.
			float gapOpen = kGapOpen;
			if (y == 0 or (y + 1 < dimY and ry[y + 1].get_seq_id() > ry[y].get_seq_id() + 1))
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

	std::vector<cif::mm::monomer *> ra, rb;

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
				if (rx[x].get_compound_id() == ry[y].get_compound_id())
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

std::tuple<std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>>
FindAtomsNearLigand(const std::vector<cif::mm::monomer *> &pa, const std::vector<cif::mm::monomer *> &pb,
	const cif::mm::residue &ra, const cif::mm::residue &rb, float maxDistance, const Ligand &ligand)
{
	float maxDistanceSq = maxDistance * maxDistance;

	std::vector<cif::mm::atom> aL, bL, aP, bP;

	std::vector<std::tuple<int, cif::mm::atom>> aI, bI;

	for (auto atom : ra.atoms())
	{
		aL.emplace_back(atom);

		int i = 0;

		for (auto &r : pa)
		{
			for (auto ra : r->atoms())
			{
				if (distance_squared(atom, ra) <= maxDistanceSq)
					aI.emplace_back(i, ra);
			}

			++i;
		}
	}

	for (auto atom : rb.atoms())
	{
		if (ligand.drops(atom.get_label_atom_id()))
			continue;

		bL.emplace_back(atom);

		int i = 0;

		for (auto &r : pb)
		{
			for (auto ra : r->atoms())
			{
				if (distance_squared(atom, ra) <= maxDistanceSq)
					bI.emplace_back(i, ra);
			}

			++i;
		}
	}

	assert(aL.size() == bL.size());
	for (size_t i = 0; i < aL.size(); ++i)
		assert(aL[i].get_label_atom_id() == ligand.map(bL[i].get_label_atom_id()));

	auto atomLess = [](const std::tuple<int, cif::mm::atom> &a, const std::tuple<int, cif::mm::atom> &b)
	{
		const auto &[ai, aa] = a;
		const auto &[bi, ba] = b;

		int d = ai - bi;

		if (d == 0)
			d = aa.get_label_atom_id().compare(ba.get_label_atom_id());

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

CAtom::CAtom(cif::atom_type type, cif::point pt, int charge, int get_seq_id, const std::string &id)
	: type(type), pt(pt), seqID(get_seq_id), id(id)
{
	const cif::atom_type_traits att(type);

	if (charge == 0)
	{
		radius = att.radius(cif::radius_type::van_der_waals);
		if (std::isnan(radius))
			radius = att.radius();
	}
	else
		radius = att.effective_ionic_radius(charge);

	if (std::isnan(radius))
		throw std::runtime_error("Unknown radius for atom " + att.symbol() + " with charge " + std::to_string(charge));
}

std::tuple<int,json> CalculateClashScore(const std::vector<CAtom> &polyAtoms, const std::vector<CAtom> &resAtoms, float maxDistance)
{
	auto maxDistanceSq = maxDistance * maxDistance;

	json distancePairs;

	int n = 0, m = 0;
	double sumOverlapSq = 0;

	for (auto &pa : polyAtoms)
	{
		bool near = false;

		for (auto &ra : resAtoms)
		{
			auto d = distance_squared(pa.pt, ra.pt);

			if (d >= maxDistanceSq)
				continue;
			
			near = true;

			d = std::sqrt(d);

			auto overlap = pa.radius + ra.radius - d;
			if (overlap < 0)
				overlap = 0;
			
			if (overlap > 0)
			{
				++n;
				sumOverlapSq += overlap * overlap;
			}
			
			json d_info{
				{ "distance", d },
				{ "VdW_overlap", overlap },
				{
					"poly_atom", {
						{ "seq_id", pa.seqID },
						{ "id", pa.id }
					}
				}
			};
			if (resAtoms.size() > 1)
				d_info["res_atom_id"] = ra.id;
			distancePairs.push_back(std::move(d_info));
		}

		if (near)
			++m;
	}

	return {
		m,
		{
			{ "score", m ? std::sqrt(sumOverlapSq / distancePairs.size()) : 0 },
			{ "clash_count", n },
			{ "poly_atom_count", m },
			{ "ligand_atom_count", resAtoms.size() },
			{ "distances", std::move(distancePairs) }
		}
	};
}

float ClashScore(cif::datablock &db, float maxDistance)
{
	using namespace cif::literals;

	auto addAtom = [](const std::string &symbol, const std::string &comp_id, int charge, cif::point p, std::vector<CAtom> &v)
	{
		const cif::atom_type_traits att(symbol);

		if (charge == 0 and att.is_metal())
		{
			auto compound = cif::compound_factory::instance().create(comp_id);
			if (compound)
				charge = compound->formal_charge();
		}

		v.emplace_back(att.type(), p, charge, 0, "");
	};

	auto &atom_site = db["atom_site"];

	std::vector<CAtom> cP, cL;

	for (const auto &[asym_id, px, py, pz, symbol, comp_id, charge] : atom_site.rows<std::string,float,float,float,std::string,std::string,int>(
			"label_asym_id", "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "label_comp_id", "pdbx_formal_charge"))
	{
		if (asym_id == "A")
			addAtom(symbol, comp_id, charge, { px, py, pz }, cP);
		else
			addAtom(symbol, comp_id, charge, { px, py, pz }, cL);
	}

	// --------------------------------------------------------------------

	auto maxDistanceSq = maxDistance * maxDistance;

	int n = 0, m = 0, o = 0;
	double sumOverlapSq = 0;

	for (auto &pa : cP)
	{
		bool near = false;

		for (auto &ra : cL)
		{
			auto d = distance_squared(pa.pt, ra.pt);

			if (d >= maxDistanceSq)
				continue;
			
			near = true;
			++o;

			d = std::sqrt(d);

			auto overlap = pa.radius + ra.radius - d;
			if (overlap < 0)
				overlap = 0;
			
			if (overlap > 0)
			{
				++n;
				sumOverlapSq += overlap * overlap;
			}
		}

		if (near)
			++m;
	}

	return m ? std::sqrt(sumOverlapSq / o) : 0;
}

// --------------------------------------------------------------------

std::tuple<float,float,float> validateCif(cif::datablock &db, const std::string &asymID, zeep::json::element &info,
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

	cif::file pdbFile(file_locator::get_pdb_file(pdbID));
	cif::mm::structure pdbStructure(pdbFile);
	cif::mm::structure afStructure(db);

	auto &ligands = LigandsTable::instance();

	for (auto hit : info["hits"])
	{
		if (hit["pdb_id"] != pdbID)
			continue;

		for (auto transplant : hit["transplants"])
		{
			if (transplant["asym_id"].as<std::string>() != asymID)
				continue;

			auto pdbAsymID = hit["pdb_asym_id"].as<std::string>();
			auto pdbget_compound_id = transplant["compound_id"].as<std::string>();

			auto ligand = ligands[pdbget_compound_id];

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

			return { rmsd1, rmsd2, rmsd3 };
		}
	}

	return { std::nanf("0"), std::nanf("0"), std::nanf("0") };
}
