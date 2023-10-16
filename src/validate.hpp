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

#include <zeep/json/element.hpp>

#include <filesystem>
#include <tuple>
#include <vector>

std::tuple<std::vector<cif::mm::monomer *>, std::vector<cif::mm::monomer *>> AlignAndTrimSequences(cif::mm::polymer &rx, cif::mm::polymer &ry);
std::tuple<std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>, std::vector<cif::mm::atom>>
FindAtomsNearLigand(const std::vector<cif::mm::monomer *> &pa, const std::vector<cif::mm::monomer *> &pb,
	const cif::mm::residue &ra, const cif::mm::residue &rb, float maxDistance, const Ligand &ligand);
std::tuple<std::vector<cif::mm::monomer *>, std::vector<cif::mm::monomer *>> AlignAndTrimSequences(cif::mm::polymer &rx, cif::mm::polymer &ry);
double Align(const cif::mm::structure &a, cif::mm::structure &b,
	std::vector<cif::point> &cAlphaA, std::vector<cif::point> &cAlphaB);
double Align(std::vector<cif::mm::atom> &aA, std::vector<cif::mm::atom> &aB);
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

std::tuple<int,zeep::json::element> CalculateClashScore(const std::vector<CAtom> &polyAtoms, const std::vector<CAtom> &resAtoms, float maxDistance);

// --------------------------------------------------------------------

float ClashScore(cif::datablock &db, float maxDistance = 4);

/// --------------------------------------------------------------------
int validateFastA(std::filesystem::path fasta, std::filesystem::path dbDir, int threads);
