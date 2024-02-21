/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2022 Maarten L. Hekkelman, NKI-AVL
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

#include "ligands.hpp"

// --------------------------------------------------------------------

std::unique_ptr<LigandsTable> LigandsTable::sInstance;

// --------------------------------------------------------------------

void Ligand::modify(cif::mm::structure &structure, const std::string &asymID) const
{
	assert(mLigand);
	auto analogue = mLigand->front()["analogue_id"].as<std::string>();
	if (not analogue.empty())
	{
		auto &res = structure.get_residue(asymID);

		std::vector<std::tuple<std::string, std::string>> remap;

		if (mModifications != nullptr)
		{
			for (const auto &[a1, a2] : mModifications->rows<std::string, std::string>("atom1", "atom2"))
				remap.emplace_back(a1, a2);
		}

		structure.change_residue(res, analogue, remap);
	}
}

size_t Ligand::atom_count(const cif::mm::residue &res) const
{
	using namespace cif::literals;

	size_t result = res.atoms().size();

	if (mModifications != nullptr)
		result -= mModifications->find("atom2"_key == cif::null).size();

	return result;
}

bool Ligand::drops(const std::string &atomID) const
{
	using namespace cif::literals;

	bool result = false;

	if (mModifications != nullptr)
		result = mModifications->contains("atom1"_key == atomID and "atom2"_key == cif::null);

	return result;
}

std::string Ligand::map(const std::string &atomID) const
{
	using namespace cif::literals;

	std::string result = atomID;

	if (mModifications != nullptr)
	{
		for (const auto &atom2 : mModifications->find<std::string>("atom1"_key == atomID, "atom2"))
		{
			result = atom2;
			break;
		}
	}
	
	return result;
}

bool LigandsTable::contains_any(const std::vector<std::string_view> &compounds) const
{
	bool result = false;

	for (auto compound_id : compounds)
	{
		Ligand l{ &mCifFile[compound_id] };
		if (l)
		{
			result = true;
			break;
		}
	}

	return result;
}
