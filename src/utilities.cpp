/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 Maarten L. Hekkelman, NKI-AVL
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

#include <filesystem>
#include <fstream>
#include <iostream>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include "revision.hpp"
#include "utilities.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

file_locator::file_locator(mcfp::config &config)
	: m_db_dir(config.get<std::string>("db-dir"))
	, m_pdb_dir(config.get<std::string>("pdb-dir"))
	, m_custom_dir(fs::path(config.get<std::string>("custom-dir")) / "out")
	, m_structure_name_pattern(config.get<std::string>("structure-name-pattern"))
	, m_pdb_name_pattern(config.get<std::string>("pdb-name-pattern"))
	, m_metadata_name_pattern(config.get<std::string>("metadata-name-pattern"))
{
	// if (not fs::is_directory(m_db_dir))
	// 	throw std::runtime_error("AlphfaFill data directory does not exist");
	// if (not fs::is_directory(m_pdb_dir))
	// 	throw std::runtime_error("PDB directory does not exist");
}

std::vector<std::filesystem::path> file_locator::get_all_structure_files(const std::string &id, int version)
{
	std::vector<fs::path> result;

	int i = 1;
	for (;;)
	{
		fs::path chunk = instance().get_structure_file(id, i, version);
		if (not fs::exists(chunk))
			break;
		
		result.emplace_back(std::move(chunk));
		++i;
	}
	
	return result;
}

// --------------------------------------------------------------------

fs::path pdbFileForID(const fs::path &pdbDir, std::string pdb_id)
{
	for (auto &ch : pdb_id)
		ch = std::tolower(ch);

	// try a PDB-REDO layout first
	fs::path pdb_path = pdbDir / pdb_id.substr(1, 2) / pdb_id / (pdb_id + "_final.cif");
	if (not fs::exists(pdb_path))
		pdb_path = pdbDir / pdb_id.substr(1, 2) / (pdb_id + ".cif.gz");

	if (not fs::exists(pdb_path))
		throw std::runtime_error("PDB file for " + pdb_id + " not found");

	return pdb_path;
}

std::vector<cif::mm::residue *> get_residuesForChain(cif::mm::structure &structure, const std::string &chain_id)
{
	std::vector<cif::mm::residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.get_asym_id() != chain_id)
			continue;

		for (auto &res : poly)
			result.emplace_back(&res);
	}

	return result;
}

using cif::point;

std::tuple<std::vector<point>, std::vector<point>> selectAtomsNearResidue(
	const std::vector<cif::mm::residue *> &pdb, const std::vector<size_t> &pdb_ix,
	const std::vector<cif::mm::residue *> &af, const std::vector<size_t> &af_ix,
	const std::vector<cif::mm::atom> &residue, float maxDistance, const Ligand &ligand)
{
	std::vector<point> ra, rb;

	float maxDistanceSq = maxDistance * maxDistance;

	assert(pdb_ix.size() == af_ix.size());

	for (size_t i = 0; i < pdb_ix.size(); ++i)
	{
		bool nearby = false;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(pdb_ix[i] < pdb.size());

			auto atom = pdb[pdb_ix[i]]->get_atom_by_atom_id(atom_id);
			if (not atom)
				continue;

			for (auto &b : residue)
			{
				if (ligand.drops(b.get_label_atom_id()))
					continue;

				if (distance_squared(atom, b) > maxDistanceSq)
					continue;

				nearby = true;
				break;
			}

			if (nearby)
				break;
		}

		if (not nearby)
			continue;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(af_ix[i] < af.size());
			assert(pdb_ix[i] < pdb.size());

			auto pt_a = pdb[pdb_ix[i]]->get_atom_by_atom_id(atom_id);
			auto pt_b = af[af_ix[i]]->get_atom_by_atom_id(atom_id);

			if (not (pt_a and pt_b))
				continue;

			ra.push_back(pt_a.get_location());
			rb.push_back(pt_b.get_location());
		}
	}

	return {ra, rb};
}


// --------------------------------------------------------------------

sequence getSequenceForStrand(cif::datablock &db, const std::string &strand)
{
	using namespace cif::literals;

	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];
	auto r = pdbx_poly_seq_scheme.find("pdb_strand_id"_key == strand);
	if (r.empty())
		throw std::runtime_error("Could not locate sequence in PDB for strand id " + strand);

	auto entity_id = r.front()["entity_id"].as<std::string>();

	auto &entity_poly = db["entity_poly"];
	auto pdb_seq = entity_poly.find1<std::string>("entity_id"_key == entity_id, "pdbx_seq_one_letter_code_can");

	return encode(pdb_seq);
}
