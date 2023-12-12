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

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "revision.hpp"
#include "utilities.hpp"

#if defined(BUILD_WEB_APPLICATION)
#include "data-service.hpp"
#endif

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::regex kAF_ID_Rx(R"((?:(AF|CS)-)?(.+?)(?:-F(\d+)(?:-(?:model|filled)_v(\d))?)?)");

std::tuple<EntryType, std::string, int, int> parse_af_id(std::string af_id)
{
#if defined(BUILD_WEB_APPLICATION)
	auto &data_service = data_service::instance();
#endif

	EntryType type = EntryType::Unknown;
	int chunkNr = 1, version = 4;
	std::string id;

	std::smatch m;
	if (std::regex_match(af_id, m, kAF_ID_Rx))
	{
		id = m[2];

		if (m[3].matched)
			chunkNr = std::stoi(m[3]);

		if (m[4].matched)
			version = std::stoi(m[4]);

		if (m[1].matched)
		{
			if (m[1] == "CS")
				type = EntryType::Custom;
			else if (m[1] == "AF")
				type = EntryType::AlphaFold;
		}
		else
		{
			// No prefix was given, try to see if we can find this ID in our cache
			for (version = 2; version < 10; ++version)
			{
				auto test = file_locator::get_metadata_file(id, chunkNr, version);

				if (fs::exists(test))
				{
					type = EntryType::AlphaFold;
					break;
				}

				if (not fs::exists(test.parent_path().parent_path()))
					break;
			}

#if defined(BUILD_WEB_APPLICATION)
			if (type != EntryType::AlphaFold and data_service.get_status(id).status != CustomStatus::Unknown)
				type = EntryType::Custom;
#endif
		}
	}

	return { type, id, chunkNr, version };
}

// --------------------------------------------------------------------

file_locator::file_locator(mcfp::config &config)
{
	if (config.has("db-dir"))
	{
		m_db_dir = config.get("db-dir");
		if (not fs::is_directory(m_db_dir))
			throw std::runtime_error("AlphaFill " + m_db_dir.string() + " data directory does not exist");
	}

	if (config.has("pdb-dir"))
	{
		m_pdb_dir = config.get("pdb-dir");
		if (not fs::is_directory(m_pdb_dir))
			throw std::runtime_error("PDB \"" + m_pdb_dir.string() + "\" directory does not exist");
	}

	if (config.has("custom-dir"))
		m_custom_dir = fs::path(config.get("custom-dir")) / "out";

	if (config.has("structure-name-pattern"))
		m_structure_name_pattern = config.get("structure-name-pattern");

	if (config.has("pdb-name-pattern"))
		m_pdb_name_pattern = config.get("pdb-name-pattern");

	if (config.has("metadata-name-pattern"))
		m_metadata_name_pattern = config.get("metadata-name-pattern");
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
		pdb_path = pdbDir / pdb_id.substr(1, 2) / pdb_id / (pdb_id + "_final.cif.gz");
	if (not fs::exists(pdb_path))
		pdb_path = pdbDir / pdb_id.substr(1, 2) / (pdb_id + ".cif.gz");

	if (not fs::exists(pdb_path))
		throw std::runtime_error("PDB file for " + pdb_id + " not found");

	return pdb_path;
}

std::vector<cif::mm::residue *> get_residuesForAsymID(cif::mm::structure &structure, const std::string &asym_id)
{
	std::vector<cif::mm::residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.get_asym_id() != asym_id)
			continue;

		for (auto &res : poly)
			result.emplace_back(&res);
	}

	return result;
}

std::vector<std::string> get_chain_ids_for_entity_id(const cif::datablock &db, const std::string &entity_id)
{
	std::vector<std::string> result;
	for (auto chain_ids : db["entity_poly"].find<std::string>(cif::key("entity_id") == entity_id, "pdbx_strand_id"))
	{
		for (auto chain_id : cif::split(chain_ids, ", ", true))
			result.emplace_back(chain_id);
	}
	return result;
}

std::vector<cif::mm::residue *> get_residues_for_chain_id(cif::mm::structure &structure, const std::string &chain_id)
{
	std::vector<cif::mm::residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.get_auth_asym_id() != chain_id)
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

		for (const char *atom_id : { "C", "CA", "N", "O" })
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

		for (const char *atom_id : { "C", "CA", "N", "O" })
		{
			assert(af_ix[i] < af.size());
			assert(pdb_ix[i] < pdb.size());

			auto pt_a = pdb[pdb_ix[i]]->get_atom_by_atom_id(atom_id);
			auto pt_b = af[af_ix[i]]->get_atom_by_atom_id(atom_id);

			if (not(pt_a and pt_b))
				continue;

			ra.push_back(pt_a.get_location());
			rb.push_back(pt_b.get_location());
		}
	}

	return { ra, rb };
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

std::vector<cif::matrix<uint8_t>> load_pae_from_file(const std::filesystem::path &file)
{
	zeep::json::element data;

	cif::gzio::ifstream in(file);

	if (not in.is_open())
		throw std::runtime_error("Could not open PAE file " + file.string());

	zeep::json::parse_json(in, data);

	if (not data.is_array())
		throw std::runtime_error("Unexpected JSON result for PAE");

	std::vector<cif::matrix<uint8_t>> result;

	for (auto e : data)
	{
		if (not e.is_object())
			throw std::runtime_error("Unexpected JSON result for PAE");
		
		auto &pae = e["predicted_aligned_error"];
		if (not pae.is_array())
			throw std::runtime_error("Unexpected JSON result for PAE");

		size_t len = pae.size();

		cif::matrix<uint8_t> &m = result.emplace_back(len, len);

		for (size_t i = 0; i < len; ++i)
		{
			for (size_t j = 0; j < len; ++j)
				m(i, j) = pae[i][j].as<int>();
		}
	}

	return result;
}
