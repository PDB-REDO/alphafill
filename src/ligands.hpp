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

#pragma once

#include <cif++/Structure.hpp>

class Ligand
{
  public:
	Ligand(cif::Datablock *db)
		: mDb(db)
		, mLigand(db ? db->get("ligand") : nullptr)
		, mModifications(db ? db->get("modification") : nullptr)
	{
	}

	Ligand(const Ligand &) noexcept = default;

	explicit operator bool() const
	{
		return mDb != nullptr and mLigand != nullptr and mLigand->front()["priority"].as<std::string>() != "n";
	}

	void modify(mmcif::Structure &structure, const std::string &asymID) const;

	std::string analogueID() const
	{
		return mLigand->front()["analogue_id"].as<std::string>();
	}

	std::string ID() const
	{
		return mDb->getName();
	}

	std::string description() const
	{
		return mLigand->front()["description"].as<std::string>();
	}

	size_t atom_count(const mmcif::Residue &res) const;

	bool drops(const std::string &atomID) const;

	std::string map(const std::string &atomID) const;

  private:
	cif::Datablock *mDb;
	cif::Category *mLigand, *mModifications;
};

// --------------------------------------------------------------------

class LigandsTable
{
  public:
	LigandsTable(const std::filesystem::path &file)
		: mCifFile(file)
	{
	}

	Ligand operator[](std::string_view id) const
	{
		return {mCifFile.get(id)};
	}

  private:
	cif::File mCifFile;
};