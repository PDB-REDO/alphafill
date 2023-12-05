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

#include "blast.hpp"
#include "ligands.hpp"

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <filesystem>
#include <regex>

// --------------------------------------------------------------------

enum class EntryType { Unknown, AlphaFold, Custom };

/// \brief Return the UniprotID and chunk number for an AlphaFold ID.
///
/// Split an id in the form of AF-UNIPROTID-F<CHUNKNR>-(model|filled)_v<VERSION>
std::tuple<EntryType,std::string,int,int> parse_af_id(std::string af_id);

// --------------------------------------------------------------------

std::filesystem::path pdbFileForID(const std::filesystem::path &pdbDir, std::string pdb_id);
std::vector<cif::mm::residue *> get_residuesForAsymID(cif::mm::structure &structure, const std::string &asym_id);
std::vector<cif::mm::residue *> get_residues_for_chain_id(cif::mm::structure &structure, const std::string &chain_id);
std::vector<std::string> get_chain_ids_for_entity_id(const cif::datablock &db, const std::string &entity_id);
std::tuple<std::vector<cif::point>, std::vector<cif::point>> selectAtomsNearResidue(
	const std::vector<cif::mm::residue *> &pdb, const std::vector<size_t> &pdb_ix,
	const std::vector<cif::mm::residue *> &af, const std::vector<size_t> &af_ix,
	const std::vector<cif::mm::atom> &residue, float maxDistance, const Ligand &ligand);
sequence getSequenceForStrand(cif::datablock &db, const std::string &strand);

extern std::regex kAF_ID_Rx;

std::vector<cif::matrix<uint8_t>> load_pae_from_file(const std::filesystem::path &file);

// --------------------------------------------------------------------

class file_locator
{
  public:
	static std::filesystem::path get_structure_file(EntryType type, const std::string &id, int chunk_nr, int version)
	{
		if (type == EntryType::Custom)
			return instance().m_custom_dir / ("CS-" + id + ".cif.gz");
		else
			return instance().get_structure_file_1(id, chunk_nr, version);
	}

	static std::filesystem::path get_metadata_file(EntryType type, const std::string &id, int chunk_nr, int version)
	{
		if (type == EntryType::Custom)
			return instance().m_custom_dir / ("CS-" + id + ".json");
		else
			return instance().get_metadata_file_1(id, chunk_nr, version);
	}

	static std::filesystem::path get_error_file(EntryType type, const std::string &id, int chunk_nr, int version)
	{
		if (type == EntryType::Custom)
			return instance().m_custom_dir / ("CS-" + id + ".json");
		else
			return instance().m_custom_dir / ("AF-" + id + "-F" + std::to_string(chunk_nr) + "-filled_v" + std::to_string(version) + ".json");
	}

	static std::filesystem::path get_structure_file(const std::string &id, int chunk_nr, int version)
	{
		return instance().get_structure_file_1(id, chunk_nr, version);
	}

	static std::filesystem::path get_pdb_file(const std::string &id)
	{
		return instance().get_pdb_file_1(id);
	}

	static std::filesystem::path get_metadata_file(const std::string &id, int chunk_nr, int version)
	{
		return instance().get_metadata_file_1(id, chunk_nr, version);
	}

	static std::vector<std::filesystem::path> get_all_structure_files(const std::string &id, int version);

  private:

	static file_locator &instance()
	{
		static file_locator s_instance(mcfp::config::instance());
		return s_instance;
	}

	file_locator(mcfp::config &config);
	file_locator(const file_locator &) = delete;
	file_locator &operator=(const file_locator &) = delete;

	std::filesystem::path get_structure_file_1(const std::string &id, int chunk_nr, int version)
	{
		std::string s = get_file(id, chunk_nr, m_structure_name_pattern);

		std::string::size_type i;
		while ((i = s.find("${db-dir}")) != std::string::npos)
			s.replace(i, strlen("${db-dir}"), m_db_dir.string());

		while ((i = s.find("${version}")) != std::string::npos)
			s.replace(i, strlen("${version}"), std::to_string(version));

		return s;
	}

	std::filesystem::path get_metadata_file_1(const std::string &id, int chunk_nr, int version)
	{
		std::string s = get_file(id, chunk_nr, m_metadata_name_pattern);

		std::string::size_type i;
		while ((i = s.find("${db-dir}")) != std::string::npos)
			s.replace(i, strlen("${db-dir}"), m_db_dir.string());
		
		while ((i = s.find("${version}")) != std::string::npos)
			s.replace(i, strlen("${version}"), std::to_string(version));

		return s;
	}

	std::filesystem::path get_pdb_file_1(const std::string &id)
	{
		std::string s = get_file(cif::to_lower_copy(id), 0, m_pdb_name_pattern);

		std::string::size_type i;
		while ((i = s.find("${pdb-dir}")) != std::string::npos)
			s.replace(i, strlen("${pdb-dir}"), m_pdb_dir.string());
		
		return s;
	}

	std::string get_file(const std::string &id, int chunk_nr, std::string pattern)
	{
		std::string::size_type i;

		std::regex rx(R"(\$\{id:(\d+):(\d+)\})");
		std::smatch m;

		while (std::regex_search(pattern, m, rx))
		{
			int start = stoi(m[1]);
			int length = stoi(m[2]);

			pattern.replace(m[0].first, m[0].second, id.substr(start, length));
		}

		while ((i = pattern.find("${id}")) != std::string::npos)
			pattern.replace(i, strlen("${id}"), id);

		while ((i = pattern.find("${chunk}")) != std::string::npos)
			pattern.replace(i, strlen("${chunk}"), std::to_string(chunk_nr));

		return pattern;
	}

	std::filesystem::path m_db_dir, m_pdb_dir, m_custom_dir;
	std::string m_structure_name_pattern;
	std::string m_pdb_name_pattern;
	std::string m_metadata_name_pattern;
};
