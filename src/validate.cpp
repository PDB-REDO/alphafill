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

#include "utilities.hpp"
#include "data-service.hpp"
#include "ligands.hpp"
#include "queue.hpp"
#include "validate.hpp"

#include <cif++.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include <fstream>
#include <iomanip>
#include <thread>

#include <cassert>

#ifdef near
#undef near
#endif

namespace fs = std::filesystem;

using json = zeep::json::element;

// --------------------------------------------------------------------

using cif::point;
using cif::quaternion;

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

double Align(const cif::mm::structure &a, cif::mm::structure &b,
	std::vector<point> &cAlphaA, std::vector<point> &cAlphaB)
{
	auto ta = center_points(cAlphaA);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << '\n';

	auto tb = center_points(cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << '\n';

	auto rotation = align_points(cAlphaB, cAlphaA);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = cif::quaternion_to_angle_axis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << '\n';
	}

	b.translate_rotate_and_translate(-tb, rotation, ta);

	for (auto &pt : cAlphaB)
		pt.rotate(rotation);

	double result = CalculateRMSD(cAlphaA, cAlphaB);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << '\n';

	return result;
}

// --------------------------------------------------------------------

std::tuple<std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>>
FindAtomsNearLigand(const std::vector<cif::mm::monomer *> &pa, const std::vector<cif::mm::monomer *> &pb,
	const cif::mm::residue &ra, const cif::mm::residue &rb, float maxDistance, const Ligand &ligand)
{
	float maxDistanceSq = maxDistance * maxDistance;

	std::vector<cif::mm::atom> aL, bL, aP, bP;

	std::vector<std::tuple<int, std::string, cif::mm::atom>> aI, bI;

	for (auto atom : ra.atoms())
	{
		aL.emplace_back(atom);

		int i = 0;

		for (auto &r : pa)
		{
			for (auto rAtom : r->atoms())
			{
				if (distance_squared(atom, rAtom) <= maxDistanceSq)
					aI.emplace_back(i, rAtom.get_label_atom_id(), rAtom);
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
			for (auto rAtom : r->atoms())
			{
				if (distance_squared(atom, rAtom) <= maxDistanceSq)
					bI.emplace_back(i, ligand.map(rAtom.get_label_atom_id()), rAtom);
			}

			++i;
		}
	}

	auto atomLess = [&ligand](const std::tuple<int, std::string, cif::mm::atom> &a, const std::tuple<int, std::string, cif::mm::atom> &b)
	{
		const auto &[ai, an, aa] = a;
		const auto &[bi, bn, ba] = b;

		int d = ai - bi;

		if (d == 0)
			d = an.compare(bn);

		return d < 0;
	};

	sort(aI.begin(), aI.end(), atomLess);
	aI.erase(std::unique(aI.begin(), aI.end()), aI.end());
	sort(bI.begin(), bI.end(), atomLess);
	bI.erase(std::unique(bI.begin(), bI.end()), bI.end());

	auto aIi = aI.begin(), bIi = bI.begin();

	while (aIi != aI.end() or bIi != bI.end())
	{
		if (aIi == aI.end())
		{
			bIi = bI.erase(bIi);
			continue;
		}

		if (bIi == bI.end())
		{
			aIi = aI.erase(aIi);
			continue;
		}

		if (atomLess(*aIi, *bIi))
		{
			aIi = aI.erase(aIi);
			continue;
		}

		if (atomLess(*bIi, *aIi))
		{
			bIi = bI.erase(bIi);
			continue;
		}

		aP.emplace_back(std::get<2>(*aIi));
		bP.emplace_back(std::get<2>(*bIi));

		++aIi;
		++bIi;
	}

	// Sort ligand atoms and make sure both are known
	std::sort(aL.begin(), aL.end(), [](const cif::mm::atom &a, const cif::mm::atom &b) { return a.get_label_atom_id().compare(b.get_label_atom_id()) < 0; });
	std::sort(bL.begin(), bL.end(), [](const cif::mm::atom &a, const cif::mm::atom &b) { return a.get_label_atom_id().compare(b.get_label_atom_id()) < 0; });

	auto aLi = aL.begin();
	auto bLi = bL.begin();

	while (aLi != aL.end() and bLi != bL.end())
	{
		auto d = aLi->get_label_atom_id().compare(bLi->get_label_atom_id());

		if (d < 0)
		{
			aLi = aL.erase(aLi);
			continue;
		}

		if (d > 0)
		{
			bLi = bL.erase(bLi);
			continue;
		}

		++aLi;
		++bLi;
	}

	if (aLi != aL.end())
		aL.erase(aLi, aL.end());
	if (bLi != bL.end())
		bL.erase(bLi, bL.end());

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

	// if (std::isnan(radius) and charge != 0)
	// {
	// 	radius = att.radius(cif::radius_type::van_der_waals);
	// 	if (std::isnan(radius))
	// 		radius = att.radius();
	// 	else
	// 		std::clog << "Error in charge for atom " << att.symbol() << " with charge " << charge << ", using charge 0 instead" << std::endl;
	// }

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
			{ "transplant_atom_count", resAtoms.size() },
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

zeep::json::element calculatePAEScore(const std::vector<cif::mm::residue *> &af_res, std::vector<CAtom> &atoms, float maxDistance, const PAE_matrix &pae)
{
	auto maxDistanceSq = maxDistance * maxDistance;

	std::vector<size_t> index;

	for (size_t i = 0; i < af_res.size(); ++i)
	{
		for (auto &a : af_res[i]->atoms())
		{
			auto l = a.get_location();
			bool near = false;

			for (auto &b : atoms)
			{
				if (distance_squared(l, b.pt) < maxDistanceSq)
				{
					near = true;
					break;
				}
			}

			if (near)
			{
				index.push_back(i);
				break;
			}
		}
	}

	zeep::json::element result;
	auto &pae_s = result["matrix"];

	std::vector<float> vt;

	for (size_t i = 0; i < index.size(); ++i)
	{
		std::vector<uint8_t> v(index.size());
		std::vector<float> vs(index.size());

		for (size_t j = 0; j < index.size(); ++j)
		{
			auto pae_v = pae(index[i], index[j]);
			
			v[j] = pae_v;

			if (i != j)
				vt.emplace_back(pae_v);
		}

		pae_s.push_back(v);
	}

	size_t N = (index.size() * (index.size() - 1));

	if (N > 1)
	{
		double sum = 0;

		for (size_t i = 0; i < index.size(); ++i)
		{
			for (size_t j = 0; j < index.size(); ++j)
			{
				if (i == j)
					continue;

				auto v = pae(index[i], index[j]);
				sum += v;
			}
		}

		double avg = sum / N;
		double sumsq = 0;

		for (size_t i = 0; i < index.size(); ++i)
		{
			for (size_t j = 0; j < index.size(); ++j)
			{
				if (i == j)
					continue;

				auto v = pae(index[i], index[j]);
				sumsq = (v - avg) * (v - avg);
			}
		}

		double stddev = std::sqrt(sumsq / N);

		result["mean"] = avg;
		result["stddev"] = stddev;

		std::sort(vt.begin(), vt.end());

		result["median"] = vt.size() % 1 == 0
			? (vt[vt.size() / 2 - 1] + vt[vt.size() / 2]) / 2.0f
			: vt[vt.size() / 2];
	}

	return result;
}

// --------------------------------------------------------------------

zeep::json::element calculateValidationScores(
	cif::datablock af_db, const std::string &asym_id,
	const std::vector<cif::mm::residue *> &pdb_res,
	const std::vector<size_t> &af_ix, const std::vector<size_t> &pdb_ix,
	const cif::mm::residue &af_ligand, const cif::mm::residue &pdb_ligand,
	float maxDistance, const Ligand &ligand)
{
	cif::mm::structure af_structure(af_db);
	auto &poly = af_structure.get_polymer_by_asym_id(asym_id);

	std::vector<cif::mm::monomer *> af_res_selected, pdb_res_selected;

	for (size_t i = 0; i < af_ix.size(); ++i)
	{
		af_res_selected.emplace_back(static_cast<cif::mm::monomer*>(&poly[af_ix[i]]));
		pdb_res_selected.emplace_back(static_cast<cif::mm::monomer*>(pdb_res[pdb_ix[i]]));

		auto &rA = *af_res_selected.back();
		auto &rP = *pdb_res_selected.back();

		if (rA.get_compound_id() != rP.get_compound_id())
			throw std::runtime_error("Invalid alignment");

		try
		{

			if (rA.get_compound_id() == "TYR" or rA.get_compound_id() == "PHE")
			{
				if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
				{
					af_structure.swap_atoms(rA.get_atom_by_atom_id("CD1"), rA.get_atom_by_atom_id("CD2"));
					af_structure.swap_atoms(rA.get_atom_by_atom_id("CE1"), rA.get_atom_by_atom_id("CE2"));
				}

				continue;
			}

			if (rA.get_compound_id() == "ASP")
			{
				if (std::abs(rA.chi(1) - rP.chi(1)) > 90)
					af_structure.swap_atoms(rA.get_atom_by_atom_id("OD1"), rA.get_atom_by_atom_id("OD2"));
				continue;
			}

			if (rA.get_compound_id() == "GLU")
			{
				if (std::abs(rA.chi(2) - rP.chi(2)) > 90)
					af_structure.swap_atoms(rA.get_atom_by_atom_id("OE1"), rA.get_atom_by_atom_id("OE2"));
				continue;
			}
		}
		catch (const std::exception &ex)
		{
			if (cif::VERBOSE > 0)
				std::cerr << ex.what() << '\n';
		}
	}

	auto &&[pA, pP, lA, lP] = FindAtomsNearLigand(af_res_selected, pdb_res_selected, af_ligand, pdb_ligand, maxDistance, ligand);

	if (pA.empty())
		return {};

	std::vector<cif::mm::atom> cA = pA, cP = pP;
	cA.insert(cA.end(), lA.begin(), lA.end());
	cP.insert(cP.end(), lP.begin(), lP.end());

	return {
		{ "local_environment_rmsd", CalculateRMSD(cA, cP) },
		{ "binding_site_rmsd", CalculateRMSD(pA, pP) },
		{ "binding_site_atom_count", pA.size() },
		{ "transplant_rmsd", CalculateRMSD(lA, lP) },
		{ "transplant_atom_count", lA.size() },
	};
}