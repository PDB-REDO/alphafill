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

#include <cif++/Cif++.hpp>
#include <cif++/CifUtils.hpp>
// #include <cif++/Compound.hpp>

#include <zeep/http/reply.hpp>
#include <zeep/json/parser.hpp>

#include "mrsrc.hpp"

#include "bsd-closefrom.h"

#include "data-service.hpp"
#include "utilities.hpp"
#include "structure.hpp"

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

#if 0
std::vector<Insertion> runBowtieInt(const std::filesystem::path& bowtie,
	const std::filesystem::path& bowtieIndex, const std::filesystem::path& fastq,
	const std::filesystem::path& logFile, unsigned threads, unsigned trimLength,
	int maxmismatch = 0, std::filesystem::path mismatchfile = {})
{
	auto p = std::to_string(threads);
	auto v = std::to_string(maxmismatch);

	std::vector<const char*> args = {
		bowtie.c_str(),
		"-m", "1",
		"-v", v.c_str(),
		"--best",
		"-p", p.c_str(),
		bowtieIndex.c_str(),
		"-"
	};

	if (maxmismatch > 0)
	{
		args.push_back("--max");
		args.push_back(mismatchfile.c_str());
	}

	args.push_back(nullptr);

	if (not fs::exists(args.front()))
		throw std::runtime_error("The executable '"s + args.front() + "' does not seem to exist");

	if (not fs::exists(fastq))
		throw std::runtime_error("The FastQ file '" + fastq.string() + "' does not seem to exist");

	// ready to roll
	int ifd[2], ofd[2], err;

	err = pipe2(ifd, O_CLOEXEC); if (err < 0) throw std::runtime_error("Pipe error: "s + strerror(errno));
	err = pipe2(ofd, O_CLOEXEC); if (err < 0) throw std::runtime_error("Pipe error: "s + strerror(errno));

	// open log file for appending
	int efd = open(logFile.c_str(), O_CREAT | O_APPEND | O_RDWR, 0644);
	const auto log_head = "\nbowtie output for " + fastq.string() + "\n" + std::string(18 + fastq.string().length(), '-') + "\n";
	write(efd, log_head.data(), log_head.size());

	int pid = fork();

	if (pid == 0)    // the child
	{
		setpgid(0, 0);        // detach from the process group, create new

		dup2(ifd[0], STDIN_FILENO);
		close(ifd[0]);
		close(ifd[1]);

		dup2(ofd[1], STDOUT_FILENO);
		close(ofd[0]);
		close(ofd[1]);

		dup2(efd, STDERR_FILENO);
		close(efd);

		closefrom(STDERR_FILENO + 1);

		const char* env[] = { nullptr };
		(void)execve(args.front(), const_cast<char* const*>(&args[0]), const_cast<char* const*>(env));
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

	std::exception_ptr ep;

	// always assume we have to trim (we used to check for trim length==read length, but that complicated the code too much)
	std::thread thread([trimLength, &fastq, fd = ifd[1], &ep]()
	{
		try
		{
			progress p(fs::file_size(fastq), fastq.string());
			p.set_action(fastq.filename().string());

			std::ifstream file(fastq, std::ios::binary);

			if (not file.is_open())
				throw std::runtime_error("Could not open file " + fastq.string());

			io::filtering_stream<io::input> in;
			std::string ext = fastq.extension().string();
			
			if (fastq.extension() == ".gz")
			{
				in.push(io::gzip_decompressor());
				ext = fastq.stem().extension().string();
			}

			in.push(progress_filter(p));
			
			in.push(file);

			char nl[1] = { '\n' };

			while (not in.eof())
			{
				// readLength != trimLength
				// read four lines

				std::string line[4];
				if (not std::getline(in, line[0]) or 
					not std::getline(in, line[1]) or 
					not std::getline(in, line[2]) or 
					not std::getline(in, line[3]))
				{
					break;
					// throw std::runtime_error("Could not read from " + fastq.string() + ", invalid file?");
				}

				if (line[0].length() < 2 or line[0][0] != '@')
					throw std::runtime_error("Invalid FastQ file " + fastq.string() + ", first line not valid");

				if (line[2].empty() or line[2][0] != '+')
					throw std::runtime_error("Invalid FastQ file " + fastq.string() + ", third line not valid");

				if (line[1].length() != line[3].length() or line[1].empty())
					throw std::runtime_error("Invalid FastQ file " + fastq.string() + ", no valid sequence data");			

				iovec v[8] = {
					{ line[0].data(), line[0].length() },
					{ nl, 1 },
					{ line[1].data(), trimLength },
					{ nl, 1 },
					{ line[2].data(), line[2].length() },
					{ nl, 1 },
					{ line[3].data(), trimLength },
					{ nl, 1 },
				};

				int r = writev(fd, v, 8);
				if (r < 0)
				{
					std::cerr << "Error writing to bowtie: " << strerror(errno) << std::endl;
					break;
				}
			}

			close(fd);
		}
		catch (const std::exception& ex)
		{
			ep = std::current_exception();
		}
	});

	close(ofd[1]);

	char buffer[8192];
	std::string line;
	std::vector<Insertion> result;

	for (;;)
	{
		int r = read(ofd[0], buffer, sizeof(buffer));

		if (r <= 0)	// keep it simple
			break;

		for (char* s = buffer; s < buffer + r; ++s)
		{
			char ch = *s;
			if (ch != '\n')
			{
				line += ch;
				continue;
			}
			
			try
			{
				auto ins = parseLine(line.c_str(), trimLength);
				if (ins.chr != INVALID)
				{
					result.push_back(ins);
					std::push_heap(result.begin(), result.end());
				}
			}
			catch (const std::exception& e)
			{
				std::cerr << std::endl
						  << "Exception parsing " << fastq << e.what() << std::endl
						  << line << std::endl
						  << std::endl;
			}

			line.clear();
		}
	}

	// should not happen... bowtie output is always terminated with a newline, right?
	if (not line.empty())
	{
		try
		{
			auto ins = parseLine(line.c_str(), trimLength);
			if (ins.chr != INVALID)
			{
				result.push_back(ins);
				std::push_heap(result.begin(), result.end());
			}
		}
		catch (const std::exception& e)
		{
			std::cerr << e.what() << std::endl
					  << line << std::endl;
		}
	}

	// return sorted and unique array of hits
	std::sort_heap(result.begin(), result.end());

	result.erase(std::unique(result.begin(), result.end()), result.end());

	thread.join();

	close(ofd[0]);
	close(efd);

	// no zombies please, removed the WNOHANG. the forked application should really stop here.
	int status = 0;
	int r = waitpid(pid, &status, 0);

	if (r == pid and WIFEXITED(status))
		status = WEXITSTATUS(status);

	if (status != 0)
		throw std::runtime_error("Error executing bowtie, result is " + std::to_string(status));

	if (ep)
		std::rethrow_exception(ep);

	return result;
}
#endif

void mergeYasaraOutput(const std::filesystem::path &input, const std::filesystem::path &yasara_out, std::ostream &os)
{
	using namespace cif::literals;

	cif::File fin(input);
	cif::File yin(yasara_out);

	auto &db_i = fin.front();
	auto &db_y = yin.front();

	auto &as_i = db_i["atom_site"];
	auto &as_y = db_y["atom_site"];

	for (auto r : as_i)
	{
		const auto &[asym_id, seq_id, atom_id, auth_seq_id] = r.get<std::string,int,std::string,int>({"label_asym_id", "label_seq_id", "label_atom_id", "auth_seq_id"});

		const auto &[x, y, z] = as_y.find1<float,float,float>(
			asym_id == "A" ?
				"label_asym_id"_key == asym_id and "label_seq_id"_key == seq_id and "label_atom_id"_key == atom_id :
				"label_asym_id"_key == asym_id and "auth_seq_id"_key == auth_seq_id and "label_atom_id"_key == atom_id,
			"Cartn_x", "Cartn_y", "Cartn_z");

		r["Cartn_x"] = x;
		r["Cartn_y"] = y;
		r["Cartn_z"] = z;
	}

	fin.save(os);
}

void optimizeWithYasara(const std::string &yasara, const std::string &af_id, std::set<std::string> requestedAsyms, float identity, std::ostream &os)
{
	using namespace std::literals;

	static std::atomic<int> sYasaraRunNr = 1;

	int my_pid = getpid();

	fs::path tmpdir = fs::temp_directory_path() / "alphafill" / (std::to_string(my_pid) + "." + std::to_string(sYasaraRunNr++));
	fs::create_directories(tmpdir);

	try
	{
		std::ofstream input(tmpdir / "input.cif");
		stripCifFile(af_id, requestedAsyms, identity, input);
		input.close();

		std::string script_s = (tmpdir / "refine.mcr").c_str();
		std::ofstream script(script_s);
		mrsrc::istream r_script("refine.mcr");
		if (not r_script)
			throw std::runtime_error("Missing resource");
		script << r_script.rdbuf();
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
			throw std::runtime_error("Error executing yasara, result is " + std::to_string(status));

		mergeYasaraOutput(tmpdir / "input.cif", tmpdir / "output.cif", os);

		fs::remove_all(tmpdir);
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << std::endl;

		throw;
	}
}
