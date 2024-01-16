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

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <zeep/http/uri.hpp>
#include <zeep/json/parser.hpp>

#include "alphafill.hpp"
#include "data-service.hpp"
#include "db-connection.hpp"
#include "https-client.hpp"
#include "queue.hpp"
#include "utilities.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

data_service &data_service::instance()
{
	static data_service s_instance;
	return s_instance;
}

data_service::data_service()
{
	auto &config = mcfp::config::instance();

	fs::path dir = config.get("custom-dir");
	m_in_dir = dir / "in";
	m_out_dir = dir / "out";
	m_work_dir = dir / "work";

	if (not fs::is_directory(m_in_dir))
		fs::create_directories(m_in_dir);

	if (not fs::is_directory(m_out_dir))
		fs::create_directories(m_out_dir);

	if (not fs::is_directory(m_work_dir))
		fs::create_directories(m_work_dir);
}

void data_service::start_queue(size_t nr_of_threads)
{
	for (size_t i = 0; i < nr_of_threads; ++i)
		m_threads.emplace_back(std::bind(&data_service::run, this));

	std::lock_guard<std::mutex> lock(m_mutex);
	std::regex rx(R"(((?:AF|CS)-.+?)(?:\.cif\.gz))");

	std::vector<std::string> files_in_in;
	for (fs::directory_iterator iter(m_in_dir); iter != fs::directory_iterator(); ++iter)
	{
		std::smatch m;

		std::string name = iter->path().filename();
		if (not std::regex_match(name, m, rx))
			continue;

		files_in_in.push_back(m[1]);
	}

	if (not files_in_in.empty())
	{
		std::thread pushThread([this, files = files_in_in]()
			{
			for (auto &file : files)
				this->m_queue.push(file); });
		pushThread.detach();
	}
}

data_service::~data_service()
{
	m_queue.push("stop");

	for (auto &t : m_threads)
	{
		if (t.joinable())
			t.join();
	}
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

void process(blocking_queue<json> &q, cif::progress_bar &p)
{
	pqxx::transaction tx1(db_connection::instance());

	for (;;)
	{
		auto data = q.pop();
		if (data.empty())
			break;

		std::string id = data["id"].as<std::string>();

		p.message(id);

		const auto &[type, uniprot_id, chunk, version] = parse_af_id(id);
		bool chunked = fs::exists(file_locator::get_metadata_file(type, uniprot_id, 2, version));

		auto r = tx1.exec1(R"(INSERT INTO af_structure (name, chunked, af_version, created, af_file) VALUES()" +
						   tx1.quote(id) + "," +
						   tx1.quote(chunked) + "," +
						   tx1.quote(data["alphafill_version"].as<std::string>()) + "," +
						   tx1.quote(data["date"].as<std::string>()) + "," +
						   tx1.quote(data["file"].as<std::string>()) +
						   ") RETURNING id");

		int64_t structure_id = r[0].as<int64_t>();

		for (auto &hit : data["hits"])
		{
			r = tx1.exec1(R"(INSERT INTO af_pdb_hit (af_id, identity, length, pdb_asym_id, pdb_id, rmsd) VALUES ()" +
						  std::to_string(structure_id) + ", " +
						  std::to_string(hit["alignment"]["identity"].as<double>()) + ", " +
						  std::to_string(hit["alignment"]["length"].as<int64_t>()) + ", " +
						  tx1.quote(hit["pdb_asym_id"].as<std::string>()) + ", " +
						  tx1.quote(hit["pdb_id"].as<std::string>()) + ", " +
						  std::to_string(hit["global_rmsd"].as<double>()) +
						  ")  RETURNING id");

			int64_t hit_id = r.front().as<int64_t>();

			for (auto &transplant : hit["transplants"])
			{
				tx1.exec0(R"(INSERT INTO af_transplant (hit_id, asym_id, compound_id, analogue_id, entity_id, rmsd) VALUES ()" +
						  std::to_string(hit_id) + ", " +
						  tx1.quote(transplant["asym_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["compound_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["analogue_id"].as<std::string>()) + ", " +
						  tx1.quote(transplant["entity_id"].as<std::string>()) + ", " +
						  std::to_string(transplant["local_rmsd"].as<double>()) +
						  ")");
			}
		}

		p.consumed(1);
	}

	tx1.commit();
}

// --------------------------------------------------------------------

int data_service::rebuild(const std::string &db_user, const fs::path &db_dir)
{
	auto &config = mcfp::config::instance();

	pqxx::work tx(db_connection::instance());

	auto schema = cif::load_resource("db-schema.sql");
	if (not schema)
		throw std::runtime_error("database schema not found (looking for db-schema.sql)");

	std::ostringstream os;
	os << schema->rdbuf();

	std::string s(os.str());

	cif::replace_all(s, "$OWNER", db_user);

	tx.exec0(s);
	tx.commit();

	std::vector<fs::path> files;
	for (auto di = fs::recursive_directory_iterator(db_dir); di != fs::recursive_directory_iterator(); ++di)
	{
		if (di->path().extension() != ".json")
			continue;
		files.push_back(di->path());
	}

	cif::progress_bar progress(files.size(), "Processing");
	blocking_queue<fs::path> q1;
	blocking_queue<json> q2;
	std::exception_ptr ep;

	std::thread t([&q2, &progress, &ep]()
		{
		try
		{
			process(q2, progress);
		}
		catch (const std::exception &ex)
		{
			ep = std::current_exception();
		} });

	std::vector<std::thread> tg;
	size_t thread_count = config.get<size_t>("threads");
	if (thread_count < 1)
		thread_count = 1;

	for (size_t i = 0; i < thread_count; ++i)
	{
		tg.emplace_back([&q1, &q2, &ep]()
			{
			for (;;)
			{
				auto f = q1.pop();

				if (f.empty())
					break;

				try
				{
					std::ifstream file(f);

					zeep::json::element data;
					zeep::json::parse_json(file, data);

					q2.push(std::move(data));
				}
				catch (const std::exception &ex)
				{
					ep = std::current_exception();
				}
			}

			q1.push({}); });
	}

	for (auto &f : files)
	{
		if (ep)
			std::rethrow_exception(ep);

		q1.push(std::move(f));
	}

	q1.push({});

	for (auto &ti : tg)
		ti.join();

	q2.push({});

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

	std::string url = config.get("alphafold-3d-beacon");

	std::string::size_type i;
	while ((i = url.find("${id}")) != std::string::npos)
		url.replace(i, strlen("${id}"), id);

	auto rep = simple_request(url);

	if (rep.get_status() == zeep::http::ok)
	{
		zeep::json::element rep_j;
		zeep::json::parse_json(rep.get_content(), rep_j);

		url = rep_j["structures"][0]["summary"]["model_url"].as<std::string>();

		rep = head_request(url, { { "Accept-Encoding", "gzip" } });

		result = rep.get_status() == zeep::http::ok;
	}

	return result;
}

std::tuple<std::filesystem::path, std::string, std::string> data_service::fetch_from_afdb(const std::string &id) const
{
	auto &config = mcfp::config::instance();

	std::string url = config.get("alphafold-3d-beacon");

	std::string::size_type i;
	while ((i = url.find("${id}")) != std::string::npos)
		url.replace(i, strlen("${id}"), id);

	auto rep = simple_request(url);

	if (rep.get_status() != zeep::http::ok)
		throw std::runtime_error("The ID " + id + " was not found at AlphaFold");

	zeep::json::element rep_j;
	zeep::json::parse_json(rep.get_content(), rep_j);

	url = rep_j["structures"][0]["summary"]["model_url"].as<std::string>();

	zeep::http::uri uri(url);

	if (uri.get_path().get_segments().empty())
		throw std::runtime_error("Empy uri returned");

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
	result << in.rdbuf();

	std::ostringstream pae_data;

	// Try PAE data
	auto o = url.find("-model_v");
	if (o != std::string::npos)
	{
		url.replace(o + 1, 5, "predicted_aligned_error");

		o = url.rfind(".cif");
		if (o != std::string::npos)
			url.replace(o + 1, 3, "json");

		rep = simple_request(url, { { "Accept-Encoding", "gzip" } });

		try
		{
			if (rep.get_status() != zeep::http::ok)
				throw std::runtime_error("Could not download PAE data");

			std::string enc = rep.get_header("Content-Encoding");
			if (enc != "gzip")
				throw std::runtime_error("Unexpected content encoding from server");

			const std::string &content = rep.get_content();

			struct membuf : public std::streambuf
			{
				membuf(char *text, size_t length)
				{
					this->setg(text, text, text + length);
				}
			} buffer(const_cast<char *>(content.data()), content.length());

			cif::gzio::istream in(&buffer);
			pae_data << in.rdbuf();
		}
		catch (const std::exception &ex)
		{
			std::cerr << "Error loading PAE data from " << std::quoted(url) << ": " << ex.what() << '\n';
		}
	}

	return { uri.get_path().get_segments().back(), result.str(), pae_data.str() };
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

// recursively print exception whats:
void print_what(std::ostream &os, const std::exception &e)
{
	os << e.what() << '\n';
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		os << " >> ";
		print_what(os, nested);
	}
}

void data_service::process_queued(const std::filesystem::path &xyzin, const std::filesystem::path &paein,
	const std::filesystem::path &xyzout, const std::filesystem::path &jsonout)
{
	std::error_code ec;

	if (not fs::exists(xyzin, ec))
		throw std::runtime_error("Input file '" + xyzin.string() + "' does not exist");

	cif::file f(xyzin);
	if (f.empty())
		throw std::runtime_error("mmCIF file seems to be empty or invalid");

	m_progress = 0;

	if (cif::VERBOSE > 0)
		std::cerr << "Running ID " << m_running << '\n';

	fs::rename(xyzin, m_work_dir / xyzin.filename(), ec);
	if (ec)
		std::cerr << "Error moving input file to work dir: " << ec.message() << '\n';

	std::vector<PAE_matrix> pae_data;

	if (not paein.empty() and fs::exists(paein))
	{
		pae_data = load_pae_from_file(paein);
		fs::rename(paein, m_work_dir / paein.filename(), ec);
	}

	auto metadata = alphafill(f.front(), pae_data, data_service_progress{ m_progress });

	try
	{
		f.save(xyzout);
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error writing output file " << std::quoted(xyzout.string()) << ": " << ex.what() << '\n';
		throw;
	}

	std::ofstream metadataFile(jsonout);
	metadataFile << metadata;

	// Clean up work files
	if (fs::exists(m_work_dir / xyzin.filename(), ec))
		fs::remove(m_work_dir / xyzin.filename(), ec);

	if (not paein.empty() and fs::exists(m_work_dir / paein.filename(), ec))
		fs::remove(m_work_dir / paein.filename(), ec);
}

void data_service::run()
{
	using namespace std::literals;
	using namespace date;
	using namespace std::chrono;

	for (;;)
	{
		std::error_code ec;
		std::filesystem::path xyzin, paein;

		auto next = m_queue.pop();

		if (next == "stop")
		{
			m_queue.push(next);
			break;
		}

		xyzin = m_in_dir / (next + ".cif.gz");

		// Skip if it does not exist
		if (not fs::exists(xyzin))
		{
			std::clog << system_clock::now() << " requested ID '" << next << "' does not resolve to filename (tried: " << xyzin.string() << ")\n";
			continue;
		}

		paein = xyzin;

		if (paein.extension() == ".gz")
			paein.replace_extension();
		if (paein.extension() == ".cif")
			paein.replace_extension("pae.gz");

		const auto &[type, afId, chunkNr, version] = parse_af_id(next);
		fs::path jsonout = file_locator::get_metadata_file(type, afId, chunkNr, version);
		fs::path xyzout = file_locator::get_structure_file(type, afId, chunkNr, version);

		if (fs::exists(xyzout, ec) and fs::exists(jsonout, ec))
		{
			// results already exist. Skip this.
			fs::remove(xyzin, ec);
			if (not paein.empty())
				fs::remove(paein, ec);
			continue;
		}

		m_running = afId;

		try
		{
			std::clog << system_clock::now() << ' ' << "process: " << xyzin.string() << '\n';
			process_queued(xyzin, paein, xyzout, jsonout);
		}
		catch (const std::exception &ex)
		{
			std::clog << system_clock::now() << ' ' << ex.what() << '\n';
			std::ofstream errorFile(m_out_dir / (next + ".error"));
			print_what(errorFile, ex);
		}

		m_progress = 0;
		m_running.clear();
	}
}

status_reply data_service::get_status(const std::string &af_id) const
{
	status_reply reply;

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
	}

	fs::path jsonFile = file_locator::get_metadata_file(type, id, chunkNr, version);
	fs::path cifFile = file_locator::get_structure_file(type, id, chunkNr, version);

	// // See if this ID might have been processed already
	// if ((not fs::exists(jsonFile) or not fs::exists(cifFile)) and not m[4].matched)
	// {
	// 	jsonFile = file_locator::get_metadata_file(type, id, chunkNr, 3);
	// 	cifFile = file_locator::get_structure_file(type, id, chunkNr, 3);
	// }

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

		std::stringstream s;
		s << in.rdbuf();

		reply.message = s.str();
	}
	else
		reply.status = CustomStatus::Unknown;

	return reply;
}

void data_service::queue(const std::string &data, const std::optional<std::string> pae, const std::string &id)
{
	std::lock_guard<std::mutex> lock(m_mutex);

	if (m_queue.is_full())
		throw std::runtime_error("The server is too busy to handle your request, please try again later");

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

	if (pae.has_value())
	{
		membuf b2(const_cast<char *>(pae->data()), pae->length());
		cif::gzio::istream in(&b2);

		cif::gzio::ofstream out(m_in_dir / (id + ".pae.gz"));
		if (not out.is_open())
			throw std::runtime_error("Could not create temporary file");

		out << in.rdbuf();
		out.close();
	}

	m_queue.push(id);
}

std::string data_service::queue_af_id(const std::string &id)
{
	if (m_queue.is_full())
		throw std::runtime_error("The server is too busy to handle your request, please try again later");

	auto &&[filename, data, pae] = fetch_from_afdb(id);

	if (filename.extension() == ".cif")
		filename.replace_extension();

	std::lock_guard<std::mutex> lock(m_mutex);

	cif::gzio::ofstream out(m_in_dir / (filename.string() + ".cif.gz"));
	if (not out.is_open())
		throw std::runtime_error("Could not create temporary file");

	out << data;
	out.close();

	if (not pae.empty())
	{
		auto paefile = m_in_dir / filename;
		paefile.replace_extension().replace_extension("pae.gz");

		cif::gzio::ofstream out(m_in_dir / paefile.filename());
		if (not out.is_open())
			throw std::runtime_error("Could not create temporary PAE file");

		out << pae;
		out.close();
	}

	// create output directory, if needed.

	const auto &[type, af_id, chunk, version] = parse_af_id(filename);
	auto outfile = file_locator::get_structure_file(af_id, chunk, version);

	assert(not outfile.empty());
	if (not fs::exists(outfile.parent_path()))
		fs::create_directories(outfile.parent_path());

	m_queue.push(filename.string());

	return filename.string();
}
