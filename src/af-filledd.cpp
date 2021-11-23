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

#include <iostream>
#include <filesystem>

#include <boost/program_options.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/algorithm/string.hpp>

#include <zeep/http/daemon.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/json/parser.hpp>
#include <zeep/crypto.hpp>
#include <zeep/http/uri.hpp>

#include <cif++/Cif++.hpp>
#include <cif++/CifUtils.hpp>

#include "db-connection.hpp"
#include "data-service.hpp"

namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;
namespace zh = zeep::http;

#define PACKAGE_NAME "af-filledd"

const uint32_t kPageSize = 20;

// --------------------------------------------------------------------

class affd_html_controller : public zh::html_controller
{
  public:
	affd_html_controller(const fs::path &db_dir)
		: mDbDir(db_dir)
	{
		mount("{,index,index.html}", &affd_html_controller::welcome);
		mount("model", &affd_html_controller::model);
		mount("structures", &affd_html_controller::structures);
		mount("compounds", &affd_html_controller::compounds);
		mount("{css,scripts,fonts,images}/", &affd_html_controller::handle_file);
		mount("favicon.ico", &affd_html_controller::handle_file);

		mount("structure-table-page", &affd_html_controller::str_table);
	}

	void welcome(const zh::request& request, const zh::scope& scope, zh::reply& reply);
	void model(const zh::request& request, const zh::scope& scope, zh::reply& reply);
	void structures(const zh::request& request, const zh::scope& scope, zh::reply& reply);
	void compounds(const zh::request& request, const zh::scope& scope, zh::reply& reply);

	void str_table(const zh::request& request, const zh::scope& scope, zh::reply& reply);

	fs::path mDbDir;
};

void affd_html_controller::welcome(const zh::request& request, const zh::scope& scope, zh::reply& reply)
{
	if (request.has_parameter("id"))
	{
		zeep::http::uri uri(request.get_uri());
		std::string afId = request.get_parameter("id");
		reply = zeep::http::reply::redirect(uri.get_path().string() + "model?id=" + zeep::http::encode_url(afId));
		return;
	}

	return get_template_processor().create_reply_from_template("index", scope, reply);
}

void affd_html_controller::structures(const zh::request& request, const zh::scope& scope, zh::reply& reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	std::string compound;
	if (request.has_parameter("compound"))
		compound = request.get_parameter("compound");

	json structures;
	auto allstructures = 
		compound.empty()
			? ds.get_structures(0, kPageSize)
			: ds.get_structures_for_compound(compound, 0, kPageSize);
	to_element(structures, allstructures);
	sub.put("structures", structures);

	sub.put("structure-count", ds.count_structures());
	sub.put("page-size", kPageSize);
	sub.put("page", 1);
	
	if (not compound.empty())
		sub.put("compound", compound);

	return get_template_processor().create_reply_from_template("structures", sub, reply);
}

void affd_html_controller::compounds(const zh::request& request, const zh::scope& scope, zh::reply& reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	json compounds;
	auto allCompounds = ds.get_compounds();
	to_element(compounds, allCompounds);
	sub.put("compounds", compounds);

	return get_template_processor().create_reply_from_template("compounds", sub, reply);
}

void affd_html_controller::model(const zh::request& request, const zh::scope& scope, zh::reply& reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	if (not request.has_parameter("id"))
		throw zeep::http::not_found;

	std::string afId = request.get_parameter("id");

	std::regex rx(R"((?:AF-)?(.+?)(?:-F1(?:-model_v1)?)?)");
	std::smatch m;
	if (std::regex_match(afId, m, rx))
		afId = m[1];

	sub.put("af_id", afId);

	fs::path jsonFile = mDbDir / ("AF-" + afId + "-F1-model_v1.cif.json");
	fs::path cifFile = mDbDir / ("AF-" + afId + "-F1-model_v1.cif.gz");

	if (not fs::exists(jsonFile) or not fs::exists(cifFile))
		throw zeep::http::not_found;

	json data;

	std::ifstream is(jsonFile);
	parse_json(is, data);

	sub.put("data", data);

	// reorder data for Tassos
	
	std::map<std::string,size_t> compoundIds;

	for (auto &hit : data["hits"])
		for (auto &transplant : hit["transplants"])
			compoundIds[transplant["compound_id"].as<std::string>()] += 1;

	json byCompound;
	for (const auto &[compoundId, count] : compoundIds)
	{
		bool firstHit = true;

		std::map<std::tuple<std::string,std::string>,size_t> transplantsPerHit;

		for (auto &hit : data["hits"])
		{
			auto key = std::make_tuple(hit["pdb_id"].as<std::string>(), hit["pdb_asym_id"].as<std::string>());
			for (auto &transplant : hit["transplants"])
			{
				if (transplant["compound_id"] == compoundId)
					transplantsPerHit[key] += 1;
			}
		}

		std::tuple<std::string,std::string> lastKey;

		for (auto &hit : data["hits"])
		{
			auto key = std::make_tuple(hit["pdb_id"].as<std::string>(), hit["pdb_asym_id"].as<std::string>());
			bool firstTransplant = key != lastKey;
			lastKey = key;

			size_t transplantCount = transplantsPerHit[key];

			for (auto &transplant : hit["transplants"])
			{
				if (transplant["compound_id"] != compoundId)
					continue;
				
				byCompound.push_back({
					{ "compound_id", compoundId },
					{ "pdb_id", hit["pdb_id"] },
					{ "pdb_asym_id", hit["pdb_asym_id"] },
					{ "alignment_length", hit["alignment_length"] },
					{ "identity", hit["identity"] },
					{ "rmsd", hit["rmsd"] },
					{ "transplant-count", transplantCount },
					{ "hit-count", count },
					{ "first-hit", firstHit },
					{ "first-transplant", firstTransplant },
					{ "transplant", transplant }
				});

				firstTransplant = false;
				firstHit = false;
			}
		}
	}

	sub.put("by_compound", byCompound);

	using namespace cif::literals;

	cif::File file(cifFile);

	std::string title = file.firstDatablock()["entity"].find1<std::string>("id"_key == 1, "pdbx_description");
	sub.put("title", title);

	get_server().get_template_processor().create_reply_from_template("model", sub, reply);
}

void affd_html_controller::str_table(const zh::request& request, const zh::scope& scope, zh::reply& reply)
{
	int page = request.get_parameter("page", 0);

	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	std::string compound;
	if (request.has_parameter("compound"))
		compound = request.get_parameter("compound");

	json structures;
	auto allstructures = 
		compound.empty()
			? ds.get_structures(page, kPageSize)
			: ds.get_structures_for_compound(compound, page, kPageSize);
	to_element(structures, allstructures);
	sub.put("structures", structures);

	return get_template_processor().create_reply_from_template("structures::structure-table-fragment", sub, reply);
}

// --------------------------------------------------------------------

class affd_rest_controller : public zh::rest_controller
{
  public:
	affd_rest_controller(const fs::path &db_dir)
		: zh::rest_controller("v1")
		, mDbDir(db_dir)
	{
		map_get_request("aff/{id}", &affd_rest_controller::get_aff_structure, "id");
		map_get_request("aff/{id}/json", &affd_rest_controller::get_aff_structure_json, "id");
	}

	zh::reply get_aff_structure(const std::string &id);
	zeep::json::element get_aff_structure_json(const std::string &id);

	fs::path mDbDir;
};

zh::reply affd_rest_controller::get_aff_structure(const std::string &id)
{
	zeep::http::reply rep(zeep::http::ok, { 1, 1 });

	fs::path file = mDbDir / ("AF-" + id + "-F1-model_v1.cif.gz");

	if (not fs::exists(file))
		return zeep::http::reply(zeep::http::not_found, { 1, 1 });
	
	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		rep.set_content(new std::ifstream(file, std::ios::binary), "text/plain");
		rep.set_header("content-encoding", "gzip");
	}
	else
	{
		std::ifstream is(file, std::ios::binary);
		if (not is.is_open())
			return zeep::http::reply(zeep::http::not_found, { 1, 1 });

		io::filtering_stream<io::input> in;
		in.push(io::gzip_decompressor());
		in.push(is);

		std::ostringstream os;
		io::copy(in, os);

		rep.set_content(os.str(), "text/plain");
	}

	rep.set_header("content-disposition", "attachement; filename = \"AF-" + id + "-F1-model_v1.cif\"");

	return rep;
}

zeep::json::element affd_rest_controller::get_aff_structure_json(const std::string &id)
{
	fs::path file = mDbDir / ("AF-" + id + "-F1-model_v1.cif.json");

	if (not fs::exists(file))
		throw zeep::http::not_found;

	zeep::json::element result;

	std::ifstream is(file);
	parse_json(is, result);

	return result;
}

// --------------------------------------------------------------------

int rebuild_db(fs::path db_dir)
{
	std::vector<fs::path> files;
	for (auto di = fs::directory_iterator(db_dir); di != fs::directory_iterator(); ++di)
	{
		if (di->path().extension() != ".json")
			continue;
		files.push_back(di->path());
	}

	cif::Progress progress(files.size(), "Processing");

	for (auto &f : files)
	{
		progress.message(f.filename());

		std::ifstream file(f);

		zeep::json::element data;
		zeep::json::parse_json(file, data);

		pqxx::transaction tx1(db_connection::instance());
		auto r = tx1.exec1(R"(INSERT INTO af_structure (name, af_version, created, af_file) VALUES()" +
			tx1.quote(data["id"].as<std::string>()) + "," +
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

		progress.consumed(1);
	}

	return 0;
}

// --------------------------------------------------------------------

int main(int argc, char* const argv[])
{
	using namespace std::literals;

	po::options_description visible("vvrd options");
	visible.add_options()
		("command",		po::value<std::string>(),	"Command, one of start, stop, status or reload")
		("no-daemon,F",								"Do not fork a background process")
		("config",		po::value<std::string>(),	"Name of config file to use, default is " PACKAGE_NAME ".conf located in current of home directory")
		("help,h",									"Show help message")
		;


	po::options_description config("config-file options");
	config.add_options()
		("address",		po::value<std::string>(),	"Address to listen to")
		("port",		po::value<unsigned short>(),"Port to listen to")
		("user",		po::value<std::string>(),	"User to run as")
		("db-dir",		po::value<std::string>(),	"Directory containing the af-filled data")

		("db-dbname",	po::value<std::string>(),	"AF DB name")
		("db-user",		po::value<std::string>(),	"AF DB owner")
		("db-password",	po::value<std::string>(),	"AF DB password")
		("db-host",		po::value<std::string>(),	"AF DB host")
		("db-port",		po::value<std::string>(),	"AF DB 5432")
		;

	po::options_description hidden("hidden options");
	hidden.add_options()
		("rebuild-db",								"Rebuild the af-filled-db")
		("debug,d", po::value<int>(),				"Debug level (for even more verbose output)");

	po::options_description cmdline_options;
	cmdline_options.add(visible).add(config).add(hidden);

	po::positional_options_description p;
	p.add("command", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = PACKAGE_NAME ".conf";
	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / PACKAGE_NAME ".conf";
	
	if (vm.count("config") != 0)
	{
		configFile = vm["config"].as<std::string>();
		if (not fs::exists(configFile))
			throw std::runtime_error("Specified config file does not seem to exist");
	}
	
	if (fs::exists(configFile))
	{
		po::options_description config_options ;
		config_options.add(config).add(hidden);

		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, config_options), vm);
	}
	
	po::notify(vm);
	int result = 0;

	try
	{
		fs::path dbDir = vm["db-dir"].as<std::string>();

		if (not fs::exists(dbDir))
			throw std::runtime_error("db-dir does not exist");

		std::vector<std::string> vConn;
		for (std::string opt: { "db-host", "db-port", "db-dbname", "db-user", "db-password" })
		{
			if (vm.count(opt) == 0)
				continue;
			
			vConn.push_back(opt.substr(3) + "=" + vm[opt].as<std::string>());
		}

		db_connection::init(ba::join(vConn, " "));

		// --------------------------------------------------------------------
		
		if (vm.count("rebuild-db"))
			return rebuild_db(dbDir);

		if (vm.count("help") or vm.count("command") == 0 or vm.count("db-dir") == 0)
		{
			std::cerr << visible << std::endl
					<< config << std::endl;
			exit(vm.count("help") ? 0 : 1);
		}

		// --------------------------------------------------------------------
		
		std::string command = vm["command"].as<std::string>();
		
		zh::daemon server([&, dbDir]()
		{
			// auto sc = new zh::security_context(secret, ibs_user_service::instance());
			// sc->add_rule("/admin", { "ADMIN" });
			// sc->add_rule("/admin/**", { "ADMIN" });
			// sc->add_rule("/media,/media/**", { "ADMIN", "EDITOR" });
			// sc->add_rule("/", {});

			auto s = new zeep::http::server(/*sc*/);

			// s->add_error_handler(new ibs_db_error_handler());

#ifndef NDEBUG
			s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#else
			s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#endif
			// s->add_controller(new zh::login_controller());
			// s->add_controller(new user_admin_rest_controller());

			s->add_controller(new affd_html_controller(dbDir));
			s->add_controller(new affd_rest_controller(dbDir));

			return s;
		}, PACKAGE_NAME);

		if (command == "start")
		{
			std::string address = "127.0.0.1";
			if (vm.count("address"))
				address = vm["address"].as<std::string>();

			unsigned short port = 10342;
			if (vm.count("port"))
				port = vm["port"].as<unsigned short>();
			
			std::string user = "www-data";
			if (vm.count("user"))
				user = vm["user"].as<std::string>();
			
			std::cout << "starting server at http://" << address << ':' << port << '/' << std::endl;

			if (vm.count("no-daemon"))
				result = server.run_foreground(address, port);
			else
				result = server.start(address, port, 1, 1, user);
		}
		else if (command == "stop")
			result = server.stop();
		else if (command == "status")
			result = server.status();
		else if (command == "reload")
			result = server.reload();
		else
		{
			std::cerr << "Invalid command" << std::endl;
			result = 1;
		}
	}
	catch (const std::exception& ex)
	{
		std::cerr << "exception:" << std::endl
			 << ex.what() << std::endl;
		result = 1;
	}

	return result;
}
