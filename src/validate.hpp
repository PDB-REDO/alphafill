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

#include "cif++/Structure.hpp"

#include "ligands.hpp"

std::tuple<float,float,float> validateCif(cif::Datablock &db, const std::string &asymID, zeep::json::element &info, float maxLigandPolyAtomDistance = 6);

std::tuple<std::vector<mmcif::Monomer *>, std::vector<mmcif::Monomer *>> AlignAndTrimSequences(mmcif::Polymer &rx, mmcif::Polymer &ry);
std::tuple<std::vector<mmcif::Atom>, std::vector<mmcif::Atom>, std::vector<mmcif::Atom>, std::vector<mmcif::Atom>>
FindAtomsNearLigand(const std::vector<mmcif::Monomer *> &pa, const std::vector<mmcif::Monomer *> &pb,
	const mmcif::Residue &ra, const mmcif::Residue &rb, float maxDistance, const Ligand &ligand);
std::tuple<std::vector<mmcif::Monomer *>, std::vector<mmcif::Monomer *>> AlignAndTrimSequences(mmcif::Polymer &rx, mmcif::Polymer &ry);
double Align(const mmcif::Structure &a, mmcif::Structure &b,
	std::vector<mmcif::Point> &cAlphaA, std::vector<mmcif::Point> &cAlphaB);
double Align(std::vector<mmcif::Atom> &aA, std::vector<mmcif::Atom> &aB);
double CalculateRMSD(const std::vector<mmcif::Atom> &a, const std::vector<mmcif::Atom> &b);

// --------------------------------------------------------------------
// clash score

struct CAtom
{
	CAtom(const CAtom &) = default;
	CAtom(CAtom &&) = default;

	CAtom(mmcif::AtomType type, mmcif::Point pt, int charge, int seqID, const std::string &id);

	CAtom(const mmcif::Atom &atom)
		: CAtom(atom.type(), atom.location(), atom.charge(), atom.labelSeqID(), atom.labelAtomID())
	{
	}

	mmcif::AtomType type;
	mmcif::Point pt;
	float radius;
	int seqID;
	std::string id;
};

std::tuple<int,zeep::json::element> CalculateClashScore(const std::vector<CAtom> &polyAtoms, const std::vector<CAtom> &resAtoms, float maxDistance);

// --------------------------------------------------------------------

float ClashScore(cif::Datablock &db, float maxDistance = 4);
