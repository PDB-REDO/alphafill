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

#include <regex>

#include <fstream>
#include <regex>
#include <thread>

#include <sys/wait.h>
#include <unistd.h>

#include <mcfp/mcfp.hpp>
#include <cif++.hpp>

#include <zeep/json/parser.hpp>

#include "mrsrc.hpp"

#include "alphafill.hpp"
#include "data-service.hpp"
#include "db-connection.hpp"
#include "https-client.hpp"
#include "queue.hpp"
#include "utilities.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

std::regex kAF_ID_Rx(R"((?:(AF|CS)-)?(.+?)(?:-F(\d+)(?:-model_v(\d))?)?)");

std::tuple<EntryType, std::string, int, int> parse_af_id(std::string af_id)
{
	auto &data_service = data_service::instance();

	EntryType type = EntryType::Unknown;
	int chunkNr = 1, version = 2;
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
			if (fs::exists(file_locator::get_metadata_file(af_id, chunkNr, 2)))
			{
				type = EntryType::AlphaFold;
				version = 2;
			}
			else if (fs::exists(file_locator::get_metadata_file(af_id, chunkNr, 3)))
			{
				type = EntryType::AlphaFold;
				version = 3;
			}
			else if (data_service.get_status(af_id).status != CustomStatus::Unknown)
				type = EntryType::Custom;
		}
	}

	return { type, id, chunkNr, version };
}

// --------------------------------------------------------------------

data_service &data_service::instance()
{
	static data_service s_instance;
	return s_instance;
}

data_service::data_service()
{
	auto &config = mcfp::config::instance();

	fs::path dir = config.get<std::string>("custom-dir");
	m_in_dir = dir / "in";
	m_out_dir = dir / "out";

	if (not fs::is_directory(m_in_dir))
		fs::create_directories(m_in_dir);

	if (not fs::is_directory(m_out_dir))
		fs::create_directories(m_out_dir);

	m_thread = std::thread(std::bind(&data_service::run, this));

	std::regex rx(R"(CS-(.+?)(?:\.cif\.gz)?)");

	for (fs::directory_iterator iter(m_in_dir); iter != fs::directory_iterator(); ++iter)
	{
		std::smatch m;

		std::string name = iter->path().filename();
		if (not std::regex_match(name, m, rx))
			continue;

		m_queue.push(m[1]);
	}
}

data_service::~data_service()
{
	m_queue.push("stop");
	m_thread.join();
}

std::vector<compound> data_service::get_compounds(float min_identity) const
{
	pqxx::work tx(db_connection::instance());

	std::vector<compound> compounds;
	for (auto const &[compound_id, analogue_id, s_count, t_count] :
		tx.stream<std::string, std::string, uint32_t, uint32_t>(
			R"(
		 select t.compound_id,
		        t.analogue_id,
				count(distinct h.af_id) as s_count,
				count(distinct t.id) as t_count
		   from af_pdb_hit h
		   join af_transplant t on h.id = t.hit_id
		  where h.identity >= )" +
			std::to_string(min_identity) + R"(
		  group by t.compound_id, t.analogue_id
	      order by t.compound_id asc)"))
	{
		compounds.emplace_back(compound{ compound_id, analogue_id, s_count, t_count });
	}

	tx.commit();

	return compounds;
}

std::vector<structure> data_service::get_structures(float min_identity, uint32_t page, uint32_t pageSize) const
{
	pqxx::work tx(db_connection::instance());

	const std::regex rx(R"(AF-(.+?-F\d+))");

	std::vector<structure> structures;
	for (auto const &[structure_id, chunked, hit_count, transplant_count, distinct] :
		tx.stream<std::string, bool, uint32_t, uint32_t, uint32_t>(
			R"(select s.name as name,
				  s.chunked,
				  count(distinct h.id) as hit_count,
				  count(distinct t.id) as transplant_count,
				  count(distinct t.analogue_id) as distinct
			 from af_structure s
			 join af_pdb_hit h on s.id = h.af_id
			 join af_transplant t on t.hit_id = h.id
			where h.identity >= )" +
			std::to_string(min_identity) + R"(
			group by s.name, s.chunked
			order by hit_count desc, s.name asc
		   offset )" +
			std::to_string(page * pageSize) + R"( rows
		    fetch first )" +
			std::to_string(pageSize) + R"( rows only)"))
	{
		const auto &[type, id, chunk, version] = parse_af_id(structure_id);
		if (chunked)
			structures.emplace_back(structure{ id + "-F" + std::to_string(chunk), hit_count, transplant_count, distinct });
		else
			structures.emplace_back(structure{ id, hit_count, transplant_count, distinct });
	}

	tx.commit();

	return structures;
}

std::vector<structure> data_service::get_structures_for_compound(float min_identity, const std::string &compound, uint32_t page, uint32_t pageSize) const
{
	pqxx::work tx(db_connection::instance());

	std::vector<structure> structures;
	for (auto const &[structure_id, chunked, hit_count, transplant_count, distinct] :
		tx.stream<std::string, bool, uint32_t, uint32_t, uint32_t>(
			R"(select s.name,
				  s.chunked,
		          count(distinct h.id) as hit_count,
				  count(distinct t.id) as transplant_count,
				  count(distinct t.analogue_id) as dist_transplant_count
			 from af_structure s
			 join af_pdb_hit h on s.id = h.af_id
			 join af_transplant t on t.hit_id = h.id
			where (t.analogue_id = )" +
			tx.quote(compound) + R"(
			   or t.compound_id = )" +
			tx.quote(compound) + R"()
			  and h.identity >= )" +
			std::to_string(min_identity) + R"(
			group by s.name, s.chunked
			order by hit_count desc
		   offset )" +
			std::to_string(page * pageSize) + R"( rows
			fetch first )" +
			std::to_string(pageSize) + R"( rows only)"))
	{
		const auto &[type, id, chunk, version] = parse_af_id(structure_id);
		if (chunked)
			structures.emplace_back(structure{ id + "-F" + std::to_string(chunk), hit_count, transplant_count, distinct });
		else
			structures.emplace_back(structure{ id, hit_count, transplant_count, distinct });
	}

	tx.commit();

	return structures;
}

uint32_t data_service::count_structures(float min_identity) const
{
	pqxx::work tx(db_connection::instance());

	auto r = tx.exec1(R"(
		  select count(distinct s.id)
		    from af_structure s
		   right join af_pdb_hit h on s.id = h.af_id
		   right join af_transplant t on t.hit_id = h.id
		   where h.identity >= )" +
					  std::to_string(min_identity));

	tx.commit();

	return r.front().as<uint32_t>();
}

uint32_t data_service::count_structures(float min_identity, const std::string &compound) const
{
	pqxx::work tx(db_connection::instance());

	auto r = tx.exec1(R"(
		  select count(distinct s.id)
		    from af_structure s
		   right join af_pdb_hit h on s.id = h.af_id
		   right join af_transplant t on t.hit_id = h.id
		   where h.identity >= )" +
					  std::to_string(min_identity) + R"(
			 and (t.compound_id = )" +
					  tx.quote(compound) + " or t.analogue_id = " + tx.quote(compound) + ")");

	tx.commit();

	return r.front().as<uint32_t>();
}

// --------------------------------------------------------------------

using json = zeep::json::element;

void process(blocking_queue<json> &q, cif::Progress &p)
{
	for (;;)
	{
		auto data = q.pop();
		if (data.empty())
			break;

		std::string id = data["id"].as<std::string>();

		p.message(id);

		const auto &[type, uniprot_id, chunk, version] = parse_af_id(id);
		bool chunked = fs::exists(file_locator::get_metadata_file(type, uniprot_id, 2, version));

		pqxx::transaction tx1(db_connection::instance());
		auto r = tx1.exec1(R"(INSERT INTO af_structure (name, chunked, af_version, created, af_file) VALUES()" +
						   tx1.quote(id) + "," +
						   tx1.quote(chunked) + "," +
						   tx1.quote(data["alphafill_version"].as<std::string>()) + "," +
						   tx1.quote(data["date"].as<std::string>()) + "," +
						   tx1.quote(data["file"].as<std::string>()) +
						   ") RETURNING id");

		int64_t structure_id = r[0].as<int64_t>();

		for (auto &alignment : data["hits"])
		{
			r = tx1.exec1(R"(INSERT INTO af_pdb_hit (af_id, identity, length, pdb_asym_id, pdb_id, rmsd) VALUES ()" +
						  std::to_string(structure_id) + ", " +
						  std::to_string(alignment["identity"].as<double>()) + ", " +
						  std::to_string(alignment["alignment_length"].as<int64_t>()) + ", " +
						  tx1.quote(alignment["pdb_asym_id"].as<std::string>()) + ", " +
						  tx1.quote(alignment["pdb_id"].as<std::string>()) + ", " +
						  std::to_string(alignment["rmsd"].as<double>()) +
						  ")  RETURNING id");

			int64_t hit_id = r.front().as<int64_t>();

			for (auto &transplant : alignment["transplants"])
			{
				tx1.exec0(R"(INSERT INTO af_transplant (hit_id, asym_id, compound_id, analogue_id, entity_id, rmsd) VALUES ()" +
						  std::to_string(hit_id) + ", " +
						  tx1.quote(transplant["asym_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["compound_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["analogue_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["entity_id"].as<std::string>()) + ", " +
						  std::to_string(alignment["rmsd"].as<double>()) +
						  ")");
			}
		}

		tx1.commit();

		p.consumed(1);
	}
}

// --------------------------------------------------------------------

int data_service::rebuild(const std::string &db_user, const fs::path &db_dir)
{
	pqxx::work tx(db_connection::instance());

#if USE_RSRC
	mrsrc::rsrc schema("db-schema.sql");
	if (not schema)
		throw std::runtime_error("database schema not found, did you include the resource?");

	std::string s(schema.data(), schema.size());

	cif::replace_all(s, "$OWNER", db_user);

	tx.exec0(s);
	tx.commit();
#else
	try
	{
		auto r = tx.exec1("select count(*) from af_structure");
		tx.commit();

		if (r.front().as<uint32_t>() > 0)
			throw std::runtime_error("not empty");
	}
	catch (const std::exception &e)
	{
		std::cerr << "Not built using resources, please create and/or empty the tables manually before running --rebuild-db" << std::endl;
		exit(1);
	}
#endif

	std::vector<fs::path> files;
	for (auto di = fs::recursive_directory_iterator(db_dir); di != fs::recursive_directory_iterator(); ++di)
	{
		if (di->path().extension() != ".json")
			continue;
		files.push_back(di->path());
	}

	cif::Progress progress(files.size(), "Processing");
	blocking_queue<json> q;
	std::exception_ptr ep;

	std::thread t([&q, &progress, &ep]()
		{
		try
		{
			process(q, progress);
		}
		catch (const std::exception &ex)
		{
			ep = std::current_exception();
		} });

	for (auto &f : files)
	{
		std::ifstream file(f);

		zeep::json::element data;
		zeep::json::parse_json(file, data);

		q.push(std::move(data));
	}

	q.push({});

	t.join();

	if (ep)
		std::rethrow_exception(ep);

	return 0;
}

// --------------------------------------------------------------------

bool data_service::exists_in_afdb(const std::string &id) const
{
	bool result = false;

	auto &config = mcfp::config::instance();

	std::string url = config.get<std::string>("alphafold-3d-beacon");

	std::string::size_type i;
	while ((i = url.find("${id}")) != std::string::npos)
		url.replace(i, strlen("${id}"), id);

	auto rep = simple_request(url);

	if (rep.get_status() == zeep::http::ok)
	{
		zeep::json::element rep_j;
		zeep::json::parse_json(rep.get_content(), rep_j);

		url = rep_j["structures"][0]["model_url"].as<std::string>();

		rep = head_request(url, { { "Accept-Encoding", "gzip" } });

		result = rep.get_status() == zeep::http::ok;
	}

	return result;
}

std::string data_service::fetch_from_afdb(const std::string &id) const
{
	auto &config = mcfp::config::instance();

	std::string url = config.get<std::string>("alphafold-3d-beacon");

	std::string::size_type i;
	while ((i = url.find("${id}")) != std::string::npos)
		url.replace(i, strlen("${id}"), id);

	auto rep = simple_request(url);

	if (rep.get_status() != zeep::http::ok)
		throw std::runtime_error("The ID " + id + " was not found at AlphaFold");

	zeep::json::element rep_j;
	zeep::json::parse_json(rep.get_content(), rep_j);

	url = rep_j["structures"][0]["model_url"].as<std::string>();

	rep = simple_request(url, { { "Accept-Encoding", "gzip" } });

	if (rep.get_status() != zeep::http::ok)
		throw std::runtime_error("Error requesting alphafold structure file: " + zeep::http::get_status_description(rep.get_status()));

	std::string enc = rep.get_header("Content-Encoding");

	if (enc != "gzip")
		throw std::runtime_error("Unexpected content encoding from server: " + enc);

	const std::string &content = rep.get_content();

	struct membuf : public std::streambuf
	{
		membuf(char *text, size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(content.data()), content.length());

	cif::gzio::istream in(&buffer);

	std::ostringstream result;

	std::string line;
	while (getline(in, line))
		result << line << std::endl;

	return result.str();
}

// --------------------------------------------------------------------

struct data_service_progress : public alphafill_progress_cb
{
	data_service_progress(std::atomic<float> &progress)
		: m_progress(progress)
	{
	}

	void set_max_0(size_t in_max) override
	{
		m_max_0 = in_max;
	}

	void set_max_1(size_t in_max) override
	{
		m_max_1 = in_max;
		++m_cur_0;
		m_cur_1 = 0;
	}

	void consumed(size_t n = 1) override
	{
		++m_cur_1;

		float p0 = static_cast<float>(m_cur_0 - 1) / m_max_0;
		float p1 = static_cast<float>(m_cur_1) / m_max_1;

		m_progress = p0 + p1 / m_max_0;
	}

	void message(const std::string &msg) override
	{
	}

	std::atomic<float> &m_progress;
	size_t m_max_0 = 1, m_max_1 = 1000, m_cur_0 = 0, m_cur_1;
};

void data_service::run()
{
	using namespace std::literals;

	for (;;)
	{
		std::string next = m_queue.pop();

		if (next == "stop")
			break;

		// std::cout << "Need to process " << next << std::endl;

		// for (int i = 0; i < 100; ++i)
		// {
		// 	m_progress = i / 100.0f;
		// 	std::this_thread::sleep_for(std::chrono::milliseconds(100));
		// }

		// m_running.clear();
		// m_progress = 0;

		std::error_code ec;

		auto xyzin = m_in_dir / (next + ".cif.gz");

		const auto &[type, afId, chunkNr, version] = parse_af_id(next);
		fs::path jsonout = file_locator::get_metadata_file(type, afId, chunkNr, version);
		fs::path xyzout = file_locator::get_structure_file(type, afId, chunkNr, version);

		if (fs::exists(xyzout, ec) and fs::exists(jsonout, ec))
		{
			// results already exist. Skip this.
			fs::remove(xyzin, ec);
			continue;
		}

		try
		{
			if (not fs::exists(xyzin, ec))
				throw std::runtime_error("Input file does not exist");

			cif::file f(xyzin);
			if (f.empty())
				throw std::runtime_error("Cif file seems to be empty or invalid");

			m_running = afId;
			m_progress = 0;

			fs::remove(xyzin, ec);

			auto metadata = alphafill(f.front(), data_service_progress{ m_progress });

			f.save(xyzout);

			std::ofstream metadataFile(jsonout);
			metadataFile << metadata;

			// int pid = fork();

			// if (pid < 0)
			// 	throw std::runtime_error("Could not fork: "s + std::strerror(errno));

			// if (pid == 0)	// child
			// {
			// 	try
			// 	{
			// 		auto metadata = alphafill(f.front(), data_service_progress{ m_progress });

			// 		f.save(xyzout);

			// 		std::ofstream metadataFile(jsonout);
			// 		metadataFile << metadata;

			// 		exit(0);
			// 	}
			// 	catch (const std::exception &ex)
			// 	{
			// 		std::ofstream errorFile(m_out_dir / ("CS-" + next + ".error"));
			// 		errorFile << ex.what() << std::endl;
			// 		exit(1);
			// 	}
			// }

			// int status = 0;
			// int err = waitpid(pid, &status, 0);

			// if (err != 0)
			// 	throw std::runtime_error("Wait failed with error "s + std::strerror(errno));

			// if (WIFEXITED(status) and WEXITSTATUS(status) != 0)
			// 	throw std::runtime_error("Alphafill terminated with exit status " + std::to_string(WEXITSTATUS(status)));

			// if (WIFSIGNALED(status))
			// 	throw std::runtime_error("Alphafill terminated with signal " + std::to_string(WTERMSIG(status)));
		}

		catch (const std::exception &ex)
		{
			std::ofstream errorFile(m_out_dir / ("CS-" + next + ".error"));
			errorFile << ex.what() << std::endl;
		}

		m_running.clear();

		if (fs::exists(xyzin, ec))
			fs::remove(xyzin, ec);
	}
}

status_reply data_service::get_status(const std::string &af_id) const
{
	status_reply reply;

	EntryType type = EntryType::Unknown;
	int chunkNr = 1, version = 2;
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
	}

	fs::path jsonFile = file_locator::get_metadata_file(type, id, chunkNr, version);
	fs::path cifFile = file_locator::get_structure_file(type, id, chunkNr, version);

	// See if this ID might have been processed already
	if ((not fs::exists(jsonFile) or not fs::exists(cifFile)) and not m[4].matched)
	{
		jsonFile = file_locator::get_metadata_file(type, id, chunkNr, 3);
		cifFile = file_locator::get_structure_file(type, id, chunkNr, 3);
	}

	if (fs::exists(jsonFile) and fs::exists(cifFile))
		reply.status = CustomStatus::Finished;
	else if (m_running == id)
	{
		reply.status = CustomStatus::Running;
		reply.progress = m_progress;
	}
	else if (fs::exists(m_in_dir / (af_id + ".cif.gz")))
		reply.status = CustomStatus::Queued;
	else if (fs::exists(m_out_dir / (af_id + ".error")))
	{
		reply.status = CustomStatus::Error;

		std::ifstream in(m_out_dir / (af_id + ".error"));

		std::string line;
		std::getline(in, line);

		reply.message = line;
	}
	else
		reply.status = CustomStatus::Unknown;

	return reply;
}

void data_service::queue(const std::string &data, const std::string &id)
{
	std::lock_guard<std::mutex> lock(m_mutex);

	struct membuf : public std::streambuf
	{
		membuf(char *text, size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(data.data()), data.length());

	cif::gzio::istream in(&buffer);

	cif::file f = cif::pdb::read(in);

	if (f.empty())
		throw std::runtime_error("Invalid or empty cif file");

	cif::gzio::ofstream out(m_in_dir / (id + ".cif.gz"));
	if (not out.is_open())
		throw std::runtime_error("Could not create temporary file");

	f.save(out);
	out.close();

	m_queue.push(id);
}

void data_service::queue_af_id(const std::string &id)
{
	std::string data = fetch_from_afdb(id);

	auto outfile = file_locator::get_metadata_file(id, 1, 3);

	cif::gzio::ofstream out(m_in_dir / ("AF-" + id + "-F1-model_v3.cif.gz"));
	if (not out.is_open())
		throw std::runtime_error("Could not create temporary file");

	out << data;
	out.close();

	// create output directory, if needed.
	assert(not outfile.empty());
	if (not fs::exists(outfile.parent_path()))
		fs::create_directories(outfile.parent_path());

	m_queue.push("AF-" + id + "-F1-model_v3");
}