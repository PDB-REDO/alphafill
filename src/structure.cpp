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

#include <cif++/Cif++.hpp>
#include <cif++/CifUtils.hpp>
// #include <cif++/Compound.hpp>

#include <zeep/json/parser.hpp>
#include <zeep/http/reply.hpp>

#include "mrsrc.hpp"

#include "utilities.hpp"
#include "data-service.hpp"

namespace fs = std::filesystem;
namespace po = boost::program_options;

// --------------------------------------------------------------------

void stripCifFile(const std::string &af_id, std::set<std::string> requestedAsyms, float identity, std::ostream &os)
{
	using namespace cif::literals;

	const auto &[id, chunkNr] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(id, chunkNr);

	if (not fs::exists(file))
		throw zeep::http::not_found;

	// optionally remove asyms whose blast origin's identity is too low
	if (identity > 0)
	{
		using json = zeep::json::element;

		fs::path jsonFile = file_locator::get_metdata_file(id, chunkNr);

		if (not fs::exists(jsonFile))
			throw zeep::http::not_found;

		json data;

		std::ifstream is(jsonFile);
		parse_json(is, data);

		for (auto &hit : data["hits"])
		{
			float hi = hit["identity"].as<float>();
			if (hi >= identity * 0.01f)
				continue;

			for (auto &transplant : hit["transplants"])
				requestedAsyms.erase(transplant["asym_id"].as<std::string>());
		}
	}

	cif::File cif(file);
	auto &db = cif.firstDatablock();

	cif.loadDictionary("mmcif_pdbx_v50");

	auto &struct_asym = db["struct_asym"];
	auto &atom_site = db["atom_site"];
	auto &struct_conn = db["struct_conn"];

	std::set<std::string> existingAsyms;
	for (const auto &[asymID] : struct_asym.rows<std::string>("id"))
		existingAsyms.insert(asymID);

	std::vector<std::string> toBeRemoved;
	std::set_difference(existingAsyms.begin(), existingAsyms.end(), requestedAsyms.begin(), requestedAsyms.end(), std::back_insert_iterator(toBeRemoved));

	for (auto &asymID : toBeRemoved)
	{
		struct_asym.erase("id"_key == asymID);
		atom_site.erase("label_asym_id"_key == asymID);
		struct_conn.erase("ptnr1_label_asym_id"_key == asymID or "ptnr2_label_asym_id"_key == asymID);
	}

	cif.save(os);
}

// --------------------------------------------------------------------

void optimizeWithYasara(const std::string &af_id, std::set<std::string> requestedAsyms, float identity, std::ostream &os)
{
	static std::atomic<int> sYasaraRunNr = 1;

	try
	{
		int pid = getpid();

		fs::path tmpdir = fs::temp_directory_path() / "alphafill" / (std::to_string(pid) + "." + std::to_string(sYasaraRunNr++));

		fs::create_directories(tmpdir);

		std::ofstream input(tmpdir / "input.cif");
		stripCifFile(af_id, requestedAsyms, identity, input);
		input.close();

		std::ofstream script(tmpdir / "refile.mcr");
		mrsrc::istream r_script("refine.mcr");
		if (not r_script)
			throw std::runtime_error("Missing resource");
		script << r_script.rdbuf();
		script.close();

	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
}
