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

#pragma once

#include "ligands.hpp"
#include "alphafill.hpp"

#include <zeep/json/element.hpp>

#include <filesystem>
#include <tuple>
#include <vector>

double Align(const cif::mm::structure &a, cif::mm::structure &b,
	std::vector<cif::point> &cAlphaA, std::vector<cif::point> &cAlphaB);

double CalculateRMSD(const std::vector<cif::mm::atom> &a, const std::vector<cif::mm::atom> &b);

// --------------------------------------------------------------------
// clash score

struct CAtom
{
	CAtom(const CAtom &) = default;
	CAtom(CAtom &&) = default;

	CAtom(cif::atom_type type, cif::point pt, int charge, int seqID, const std::string &id);

	CAtom(const cif::mm::atom &atom)
		: CAtom(atom.get_type(), atom.get_location(), atom.get_charge(), atom.get_label_seq_id(), atom.get_label_atom_id())
	{
	}

	cif::atom_type type;
	cif::point pt;
	float radius;
	int seqID;
	std::string id;
};

std::tuple<int, zeep::json::element> CalculateClashScore(const std::vector<CAtom> &polyAtoms, const std::vector<CAtom> &resAtoms, float maxDistance);

float ClashScore(cif::datablock &db, float maxDistance = 4);

// --------------------------------------------------------------------

zeep::json::element calculateValidationScores(
	cif::datablock af_db, const std::string &asym_id,
	const std::vector<cif::mm::residue *> &pdb_res,
	const std::vector<size_t> &af_ix, const std::vector<size_t> &pdb_ix,
	const cif::mm::residue &af_ligand, const cif::mm::residue &pdb_ligand,
	float maxDistance, const Ligand &ligand);

// --------------------------------------------------------------------

zeep::json::element calculatePAEScore(const std::vector<cif::mm::residue *> &af_res, std::vector<CAtom> &atoms, float maxDistance, const PAE_matrix &pae);
