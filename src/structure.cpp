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

#include <fcntl.h>
#include <fstream>
#include <future>
#include <regex>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/wait.h>
#include <unistd.h>

#include <atomic>
#include <filesystem>
#include <fstream>

#include <cif++.hpp>

#include <zeep/http/reply.hpp>
#include <zeep/json/parser.hpp>

#include "bsd-closefrom.h"

#include "data-service.hpp"
#include "utilities.hpp"
#include "structure.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;

using json = zeep::json::element;

// --------------------------------------------------------------------

void stripCifFile(const std::string &af_id, std::set<std::string> requestedAsyms, float identity, std::ostream &os)
{
	using namespace cif::literals;

	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(type, id, chunkNr, version);

	if (not fs::exists(file))
		throw zeep::http::not_found;

	// optionally remove asyms whose blast origin's identity is too low
	if (identity > 0)
	{
		using json = zeep::json::element;

		fs::path jsonFile = file_locator::get_metadata_file(type, id, chunkNr, version);

		if (not fs::exists(jsonFile))
			throw zeep::http::not_found;

		json data;

		std::ifstream is(jsonFile);
		parse_json(is, data);

		for (auto &hit : data["hits"])
		{
			float hi = hit["alignment"]["identity"].as<float>();
			if (hi >= identity * 0.01f)
				continue;

			for (auto &transplant : hit["transplants"])
				requestedAsyms.erase(transplant["asym_id"].as<std::string>());
		}
	}

	cif::file cif(file);
	auto &db = cif.front();

	if (cif.get_validator() == nullptr)
		cif.load_dictionary("mmcif_af");

	auto &struct_asym = db["struct_asym"];
	auto &atom_site = db["atom_site"];
	auto &struct_conn = db["struct_conn"];
	auto &entity_poly = db["entity_poly"];

	std::set<std::string> existingAsyms;
	for (const auto &[asymID, entityID] : struct_asym.rows<std::string,std::string>("id", "entity_id"))
	{
		// check if this is a nonpoly entity
		if (entity_poly.exists("entity_id"_key == entityID))
			continue;

		existingAsyms.insert(asymID);
	}

	std::vector<std::string> toBeRemoved;
	std::set_difference(existingAsyms.begin(), existingAsyms.end(), requestedAsyms.begin(), requestedAsyms.end(), std::back_insert_iterator(toBeRemoved));

	auto validator = db.get_validator();
	db.set_validator(nullptr);

	for (auto &asymID : toBeRemoved)
	{
		struct_asym.erase("id"_key == asymID);
		atom_site.erase("label_asym_id"_key == asymID);
		struct_conn.erase("ptnr1_label_asym_id"_key == asymID or "ptnr2_label_asym_id"_key == asymID);
	}

	// fix up pdbx_struct_assembly_gen, if it exists?
	if (db.get("pdbx_struct_assembly_gen") != nullptr)
	{
		auto &pdbx_struct_assembly_gen = db["pdbx_struct_assembly_gen"];
		for (auto r : pdbx_struct_assembly_gen)
		{
			auto asym_id_list = cif::split<std::string>(r["asym_id_list"].as<std::string>(), ",", true);

			std::vector<std::string> new_asym_id_list;
			std::set_intersection(asym_id_list.begin(), asym_id_list.end(), requestedAsyms.begin(), requestedAsyms.end(), std::back_insert_iterator(new_asym_id_list));
			r["asym_id_list"] = cif::join(new_asym_id_list, ",");
		}
	}

	db.set_validator(validator);

	cif::mm::structure structure(db);
	structure.cleanup_empty_categories();

	cif.save(os);
}

// --------------------------------------------------------------------

json mergeYasaraOutput(const std::filesystem::path &input, const std::filesystem::path &yasara_out, std::ostream &os)
{
	using namespace cif::literals;

	cif::file fin(input);
	cif::file yin(yasara_out);

	auto &db_i = fin.front();
	auto &db_y = yin.front();

	json info;
	const auto &[type, afID, chunkNr, version] = parse_af_id(db_i.name());
	std::ifstream infoFile(file_locator::get_metadata_file(type, afID, chunkNr, version));
	zeep::json::parse_json(infoFile, info);

	// statistics before

	std::string ligandAsymID = db_i["struct_asym"].find1<std::string>("id"_key != "A", "id");

	float clashBefore = ClashScore(fin.front());

	auto &as_y = db_y["atom_site"];

	using key_type = std::tuple<std::string,int,std::string>;
	using value_type = std::tuple<float,float,float>;
	std::map<key_type, value_type> locations;

	for (const auto &[asym_id, seq_id, atom_id, auth_seq_id, x, y, z] :
		as_y.find<std::string,int,std::string,int,float,float,float>("type_symbol"_key != "H",
			"label_asym_id", "label_seq_id", "label_atom_id", "auth_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"))
	{
		if (asym_id == "A")
			locations.emplace(key_type{ asym_id, seq_id, atom_id }, value_type{ x, y, z });
		else
			locations.emplace(key_type{ asym_id, 0, atom_id }, value_type{ x, y, z });
	}

	auto &as_i = db_i["atom_site"];

	for (auto r : as_i.rows())
	{
		const auto &[asym_id, seq_id, atom_id, auth_seq_id] = r.get<std::string,int,std::string,int>("label_asym_id", "label_seq_id", "label_atom_id", "auth_seq_id");

		auto l = locations.find(asym_id == "A" ? key_type{ asym_id, seq_id, atom_id } : key_type{ asym_id, 0, atom_id });
		if (l == locations.end())
			continue;

		r["Cartn_x"] = std::get<0>(l->second);
		r["Cartn_y"] = std::get<1>(l->second);
		r["Cartn_z"] = std::get<2>(l->second);
	}

	fin.save(os);

	float clashAfter = ClashScore(fin.front());

	return{
		{ "clash", {
			{ "before", clashBefore },
			{ "after", clashAfter }
		}}
	};
}

json optimizeWithYasara(const std::string &af_id, std::set<std::string> requestedAsyms, std::ostream &os)
{
	using namespace std::literals;

	auto &config = mcfp::config::instance();
	auto yasara = config.get("yasara");

	static std::atomic<int> sYasaraRunNr = 1;

	int my_pid = getpid();

	fs::path tmpdir = fs::temp_directory_path() / "alphafill" / (std::to_string(my_pid) + "." + std::to_string(sYasaraRunNr++));
	fs::create_directories(tmpdir);

	std::ofstream input(tmpdir / "input.cif");
	stripCifFile(af_id, requestedAsyms, 0, input);
	input.close();

	std::string script_s = (tmpdir / "refine.mcr").c_str();
	std::ofstream script(script_s);
	auto r_script = cif::load_resource("refine.mcr");
	if (not r_script)
		throw std::runtime_error("Missing resource refine.mcr");
	script << r_script->rdbuf();
	script.close();

	std::string modelin = ("modelin='" + (tmpdir / "input.cif").string() + "'");
	std::string modelout = ("modelout='" + (tmpdir / "output.cif").string() + "'");

	std::vector<const char *> args = {
		yasara.c_str(),
		"-txt",
		script_s.c_str(),
		modelin.c_str(),
		modelout.c_str(),
		nullptr};

	if (not fs::exists(args.front()))
		throw std::runtime_error("The executable '"s + args.front() + "' does not seem to exist");

	// ready to roll
	int ifd[2], ofd[2], err;

	err = pipe2(ifd, O_CLOEXEC);
	if (err < 0)
		throw std::runtime_error("Pipe error: "s + strerror(errno));
	err = pipe2(ofd, O_CLOEXEC);
	if (err < 0)
		throw std::runtime_error("Pipe error: "s + strerror(errno));

	// open log file for appending
	int efd = open((tmpdir / "yasara.log").c_str(), O_CREAT | O_APPEND | O_RDWR, 0644);
	const auto log_head = "\nyasara output for " + af_id + "\n" + std::string(18 + af_id.length(), '-') + "\n";
	write(efd, log_head.data(), log_head.size());

	int pid = fork();

	if (pid == 0) // the child
	{
		setpgid(0, 0); // detach from the process group, create new

		dup2(ifd[0], STDIN_FILENO);
		close(ifd[0]);
		close(ifd[1]);

		dup2(ofd[1], STDOUT_FILENO);
		close(ofd[0]);
		close(ofd[1]);

		dup2(efd, STDERR_FILENO);
		close(efd);

		closefrom(STDERR_FILENO + 1);

		std::error_code ec;
		fs::current_path(tmpdir, ec);

		const char *env[] = {nullptr};
		(void)execve(args.front(), const_cast<char *const *>(&args[0]), const_cast<char *const *>(env));
		exit(-1);
	}

	if (pid == -1)
	{
		close(ifd[0]);
		close(ifd[1]);
		close(ofd[0]);
		close(ofd[1]);
		close(efd);

		throw std::runtime_error("fork failed: "s + strerror(errno));
	}

	close(ifd[0]);
	close(ofd[1]);

	// start reading output

	char buffer[8192];
	std::string line;
	bool done = false;

	for (;;)
	{
		int r = read(ofd[0], buffer, sizeof(buffer));

		if (r <= 0) // keep it simple
			break;

		write(efd, buffer, r);

		for (char *s = buffer; s < buffer + r; ++s)
		{
			char ch = *s;
			if (ch != '\n')
			{
				line += ch;
				continue;
			}

			if (line.substr(0, 11) == " - ERROR - ")
				break;

			if (line.substr(0, 4) == "DONE")
			{
				done = true;
				break;
			}

			line.clear();
		}
	}

	close(ofd[0]);
	close(efd);

	if (not done)
		kill(pid, 9);

	// no zombies please, removed the WNOHANG. the forked application should really stop here.
	int status = 0;
	int r = waitpid(pid, &status, 0);

	if (r == pid and WIFEXITED(status))
		status = WEXITSTATUS(status);

	if (status != 0)
	{
		std::stringstream msg;
		msg << "Error executing yasara, exit code is: " << status << '\n';

		if (fs::exists(tmpdir / "errorexit.txt"))
		{
			try
			{
				std::ifstream errormsg(tmpdir / "errorexit.txt");
				msg << errormsg.rdbuf();
			}
			catch (...)
			{
			}
		}

		throw std::runtime_error(msg.str());
	}

	auto score = mergeYasaraOutput(tmpdir / "input.cif", tmpdir / "output.cif", os);

	fs::remove_all(tmpdir);

	return score;
}
