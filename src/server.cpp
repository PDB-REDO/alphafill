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
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <zeep/crypto.hpp>
#include <zeep/http/daemon.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/http/uri.hpp>
#include <zeep/json/parser.hpp>

#include <cif++.hpp>
#include <cfg.hpp>

#include "data-service.hpp"
#include "db-connection.hpp"
#include "ligands.hpp"
#include "server.hpp"
#include "structure.hpp"
#include "utilities.hpp"

namespace ba = boost::algorithm;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace zh = zeep::http;

#define PACKAGE_NAME "af-filledd"

const std::array<uint32_t, 6> kIdentities{70, 60, 50, 40, 30, 25};

const int
	kPageSize = 20,
	kMinIdentity = kIdentities.back(),
	kMaxIdentity = kIdentities.front();

// --------------------------------------------------------------------
// Register an object to handle #af.link('pdb-id') calls from the template
// processor.

class af_link_template_object : public zh::expression_utility_object<af_link_template_object>
{
  public:
	static constexpr const char *name() { return "af"; }

	void set_template(const std::string &t)
	{
		m_template = t;
	}

	virtual zh::object evaluate(const zh::scope &scope, const std::string &methodName,
		const std::vector<zh::object> &parameters) const
	{
		zh::object result;

		if (methodName == "link" and parameters.size() >= 1 and not m_template.empty())
		{
			try
			{
				auto id = parameters.front().as<std::string>();
				std::regex rx(R"([0-9][0-9a-z]{3}(?:\.[a-z]+))", std::regex_constants::icase);

				if (std::regex_match(id, rx))
					id.erase(4, std::string::npos);

				auto url = m_template;
				std::string::size_type s;
				while ((s = url.find("${id}")) != std::string::npos)
					url.replace(s, 5, id);

				to_element(result, url);
			}
			catch (const std::exception &e)
			{
				std::cerr << "Error getting post by ID: " << e.what() << std::endl;
			}
		}

		return result;
	}

  private:
	std::string m_template;
} s_af_object;

// --------------------------------------------------------------------

class missing_entry_error : public std::runtime_error
{
  public:
	missing_entry_error(const std::string &id)
		: runtime_error(id)
	{
	}
};

// --------------------------------------------------------------------

class missing_entry_error_handler : public zeep::http::error_handler
{
  public:
	virtual bool create_error_reply(const zeep::http::request &req, std::exception_ptr eptr, zeep::http::reply &reply);
};

bool missing_entry_error_handler::create_error_reply(const zeep::http::request &req, std::exception_ptr eptr, zeep::http::reply &reply)
{
	bool result = false;

	try
	{
		std::rethrow_exception(eptr);
	}
	catch (const missing_entry_error &ex)
	{
		zeep::http::scope scope(*get_server(), req);
		scope.put("missing", ex.what());

		get_server()->get_template_processor().create_reply_from_template("index", scope, reply);
		result = true;
	}
	catch (...)
	{
	}

	return result;
}

// --------------------------------------------------------------------

class affd_html_controller : public zh::html_controller
{
  public:
	affd_html_controller()
	{
		mount("{,index,index.html}", &affd_html_controller::welcome);
		mount("model", &affd_html_controller::model);
		mount("optimized", &affd_html_controller::optimized);
		mount("structures", &affd_html_controller::structures);
		mount("compounds", &affd_html_controller::compounds);
		mount("about", &affd_html_controller::about);
		mount("download", &affd_html_controller::download);
		mount("{css,scripts,fonts,images}/", &affd_html_controller::handle_file);
		mount("alphafill.json.schema", &affd_html_controller::schema);
		mount("favicon.ico", &affd_html_controller::handle_file);

		mount("structure-table-page", &affd_html_controller::structures_table);
	}

	void welcome(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void model(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void optimized(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void structures(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void compounds(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void about(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void download(const zh::request &request, const zh::scope &scope, zh::reply &reply);
	void schema(const zh::request &request, const zh::scope &scope, zh::reply &reply);

	void structures_table(const zh::request &request, const zh::scope &scope, zh::reply &reply);
};

void affd_html_controller::welcome(const zh::request &request, const zh::scope &scope, zh::reply &reply)
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

void affd_html_controller::structures(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	std::string compound;
	if (request.has_parameter("compound"))
		compound = request.get_parameter("compound");

	int identity = kMaxIdentity;
	if (request.has_parameter("identity"))
		identity = std::stoi(request.get_parameter("identity"));

	if (identity < kMinIdentity)
		identity = kMinIdentity;
	if (identity > 100)
		identity = 100;

	sub.put("identity", identity);

	json structures;
	auto allstructures =
		compound.empty()
			? ds.get_structures(identity * 0.01f, 0, kPageSize)
			: ds.get_structures_for_compound(identity * 0.01f, compound, 0, kPageSize);
	to_element(structures, allstructures);
	sub.put("structures", structures);

	sub.put("structure-count",
		compound.empty()
			? ds.count_structures(identity * 0.01f)
			: ds.count_structures(identity * 0.01f, compound));
	sub.put("page-size", kPageSize);
	sub.put("page", 1);

	if (not compound.empty())
		sub.put("compound", compound);

	return get_template_processor().create_reply_from_template("structures", sub, reply);
}

void affd_html_controller::structures_table(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	int page = request.get_parameter("page", 0);

	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	std::string compound;
	if (request.has_parameter("compound"))
		compound = request.get_parameter("compound");

	int identity = kMaxIdentity;
	if (request.has_parameter("identity"))
		identity = std::stoi(request.get_parameter("identity"));

	if (identity < kMinIdentity)
		identity = kMinIdentity;
	if (identity > 100)
		identity = 100;

	sub.put("identity", identity);

	json structures;
	auto allstructures =
		compound.empty()
			? ds.get_structures(identity * 0.01f, page, kPageSize)
			: ds.get_structures_for_compound(identity * 0.01f, compound, page, kPageSize);
	to_element(structures, allstructures);
	sub.put("structures", structures);

	return get_template_processor().create_reply_from_template("structures::structure-table-fragment", sub, reply);
}

void affd_html_controller::compounds(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	auto &ds = data_service::instance();

	int identity = kMaxIdentity;
	if (request.has_parameter("identity"))
		identity = std::stoi(request.get_parameter("identity"));

	if (identity < kMinIdentity)
		identity = kMinIdentity;
	if (identity > 100)
		identity = 100;

	sub.put("identity", identity);

	json compounds;
	auto allCompounds = ds.get_compounds(identity * 0.01f);
	to_element(compounds, allCompounds);
	sub.put("compounds", compounds);

	return get_template_processor().create_reply_from_template("compounds", sub, reply);
}

struct transplant_info
{
	std::string compound_id, analogue_id;
	int hit_nr;
	std::string pdb_id;
	double identity;
	double gRMSd;
	std::string asym_id;
	double clashScore;
	double lRMSd;
	bool firstHit = false;
	bool firstTransplant = false;
	int hitCount = 1;
	int transplantCount = 1;

	template <typename Archive>
	void serialize(Archive &ar, unsigned long)
	{
		ar &zeep::make_nvp("compound_id", compound_id) & zeep::make_nvp("analogue_id", analogue_id) & zeep::make_nvp("pdb_id", pdb_id) & zeep::make_nvp("identity", identity) & zeep::make_nvp("global-rmsd", gRMSd) & zeep::make_nvp("asym_id", asym_id) & zeep::make_nvp("local-rmsd", lRMSd) & zeep::make_nvp("clash-score", clashScore) & zeep::make_nvp("first-hit", firstHit) & zeep::make_nvp("first-transplant", firstTransplant) & zeep::make_nvp("hit-count", hitCount) & zeep::make_nvp("transplant-count", transplantCount);
	}

	bool operator<(const transplant_info &rhs) const
	{
		int d = compound_id.compare(rhs.compound_id);

		if (d == 0)
		{
			double dd = gRMSd - rhs.gRMSd;
			if (dd != 0)
				d = std::signbit(dd) ? -1 : 1;
		}

		if (d == 0)
			d = pdb_id.compare(rhs.pdb_id);

		if (d == 0)
		{
			double dd = lRMSd - rhs.lRMSd;
			if (dd != 0)
				d = std::signbit(dd) ? -1 : 1;
		}

		if (d == 0)
			d = asym_id.compare(rhs.asym_id);

		return d < 0;
	}
};

void affd_html_controller::model(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	if (not request.has_parameter("id"))
		throw missing_entry_error("<missing-id>");

	const auto &[type, afId, chunkNr] = parse_af_id(request.get_parameter("id"));

	if (type == EntryType::Custom)
	{
		const auto &[status, progress] = data_service::instance().get_status(afId);
		if (status == CustomStatus::Queued or status == CustomStatus::Running)
		{
			sub.put("hash", afId);
			sub.put("status", status);
			get_server().get_template_processor().create_reply_from_template("wait", sub, reply);
			return;
		}
	}

	sub.put("af_id", afId);
	sub.put("chunk", chunkNr);

	bool chunked = chunkNr > 1 or fs::exists(file_locator::get_metadata_file(type, afId, 2));

	sub.put("chunked", chunked);

	if (chunked)
	{
		auto allChunks = file_locator::get_all_structure_files(afId);
		json chunks;

		for (size_t i = 0; i < allChunks.size(); ++i)
			chunks.emplace_back(afId + "-F" + std::to_string(i + 1));
		
		sub.put("chunks", chunks);
	}

	fs::path jsonFile = file_locator::get_metadata_file(type, afId, chunkNr);
	fs::path cifFile = file_locator::get_structure_file(type, afId, chunkNr);

	if (not fs::exists(jsonFile) /*or not fs::exists(cifFile)*/)
		throw missing_entry_error(afId);

	json data;

	std::ifstream is(jsonFile);
	parse_json(is, data);

	int identity;

	if (request.has_parameter("identity"))
	{
		identity = std::stoi(request.get_parameter("identity"));
		if (identity < kMinIdentity)
			identity = kMinIdentity;
		if (identity > 100)
			identity = 100;
	}
	else // see article, take highest identity that results in transplants
	{
		identity = kMinIdentity;
		for (int i : kIdentities)
		{
			for (auto &hit : data["hits"])
			{
				if (hit["identity"].as<double>() * 100 < i)
					continue;

				identity = i;
				break;
			}

			if (identity > kMinIdentity)
				break;
		}
	}

	sub.put("identity", identity);

	int hit_nr = 0;
	std::vector<transplant_info> transplants;
	for (auto &hit : data["hits"])
	{
		double hitIdentity = hit["identity"].as<double>();
		if (hitIdentity * 100 < identity)
			continue;

		++hit_nr;
		for (auto &transplant : hit["transplants"])
		{
			transplants.emplace_back(transplant_info{
				transplant["compound_id"].as<std::string>(),
				transplant["analogue_id"].as<std::string>(),
				hit_nr,
				hit["pdb_id"].as<std::string>() + '.' + hit["pdb_asym_id"].as<std::string>(),
				hitIdentity,
				hit["rmsd"].as<double>(),
				transplant["asym_id"].as<std::string>(),
				transplant["clash"]["score"].as<double>(),
				transplant["rmsd"].as<double>()});
		}
	}

	std::sort(transplants.begin(), transplants.end());

	auto e = transplants.end();
	for (auto i = transplants.begin(); i != transplants.end();)
	{
		i->firstHit = true;

		auto j = i + 1;
		while (j != e and j->compound_id == i->compound_id)
			++j;

		size_t n = j - i;
		for (; i != j; ++i)
			i->hitCount = n;
	}

	for (auto i = transplants.begin(); i != e;)
	{
		i->firstTransplant = true;

		auto j = i + 1;
		while (j != e and j->compound_id == i->compound_id and j->hit_nr == i->hit_nr)
			++j;

		size_t n = j - i;
		for (; i != j; ++i)
			i->transplantCount = n;
	}

	json byCompound;
	to_element(byCompound, transplants);

	sub.put("by_compound", byCompound);

	using namespace cif::literals;

	try
	{
		cif::file file(cifFile);
		if (file.empty())
			throw zeep::http::not_found;

		std::string title = file.front()["entity"].find1<std::string>("id"_key == 1, "pdbx_description");
		sub.put("title", title);
	}
	catch (const std::exception &e)
	{
		sub.put("title", e.what());
	}

	// TODO: These magic numbers should of course be configurable parameters
	// 11.43, 3.04 voor global, en 3.10 en 0.92 voor local
	sub.put("cutoff", json{
		{"global", {{"unreliable", 11.43}, {"suspect", 3.04}}},
		{"local", {{"unreliable", 3.10}, {"suspect", 0.92}}}
	});

	get_server().get_template_processor().create_reply_from_template("model", sub, reply);
}

void affd_html_controller::optimized(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	using json = zeep::json::element;

	zh::scope sub(scope);

	if (not request.has_parameter("id"))
		throw missing_entry_error("<missing-id>");

	if (not request.has_parameter("asym"))
		throw missing_entry_error("<missing-asym>");

	std::string asymID = request.get_parameter("asym");

	const auto &[type, afId, chunkNr] = parse_af_id(request.get_parameter("id"));

	sub.put("af_id", afId);
	sub.put("chunk", chunkNr);
	sub.put("asym_id", asymID);

	try
	{
		using namespace cif::literals;

		auto cif = file_locator::get_structure_file(type, afId, chunkNr);
		cif::file cf(cif);
		if (cf.empty())
			throw zeep::http::not_found;

		auto &db = cf.front();
		auto entity_id = db["struct_asym"].find1<std::string>("id"_key == asymID, "entity_id");
		const auto &[compound_name, compound_id] = db["pdbx_entity_nonpoly"].find1<std::string, std::string>("entity_id"_key == entity_id, "name", "comp_id");

		sub.put("compound-name", compound_name);
		sub.put("compound-id", compound_id);
	}
	catch (...)
	{
	}

	bool chunked = chunkNr > 1 or fs::exists(file_locator::get_metadata_file(type, afId, 2));

	sub.put("chunked", chunked);

	if (chunked)
	{
		auto allChunks = file_locator::get_all_structure_files(afId);
		json chunks;

		for (size_t i = 0; i < allChunks.size(); ++i)
			chunks.emplace_back(afId + "-F" + std::to_string(i + 1));
		
		sub.put("chunks", chunks);
	}

	get_server().get_template_processor().create_reply_from_template("optimized", sub, reply);
}

void affd_html_controller::about(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	return get_template_processor().create_reply_from_template("about", scope, reply);
}

void affd_html_controller::download(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	return get_template_processor().create_reply_from_template("download", scope, reply);
}

void affd_html_controller::schema(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	html_controller::handle_file(request, scope, reply);
	if (reply.get_status() == zeep::http::ok)
		reply.set_header("Content-Disposition", R"(attachment; filename="alphafill.json.schema")");
}

// --------------------------------------------------------------------

class affd_rest_controller : public zh::rest_controller
{
  public:
	affd_rest_controller()
		: zh::rest_controller("v1")
	{
		map_get_request("aff/{id}", &affd_rest_controller::get_aff_structure, "id");
		map_get_request("aff/{id}/json", &affd_rest_controller::get_aff_structure_json, "id");
		map_get_request("aff/{id}/status", &affd_rest_controller::get_aff_status, "id");

		map_get_request("aff/3d-beacon/{id}", &affd_rest_controller::get_aff_3d_beacon, "id");

		map_get_request("aff/{id}/stripped/{asymlist}", &affd_rest_controller::get_aff_structure_stripped_def, "id", "asymlist");
		map_get_request("aff/{id}/stripped/{asymlist}/{identity}", &affd_rest_controller::get_aff_structure_stripped, "id", "asymlist", "identity");

		map_get_request("aff/{id}/optimized/{asymlist}", &affd_rest_controller::get_aff_structure_optimized, "id", "asymlist");
		map_get_request("aff/{id}/optimized-with-stats/{asymlist}", &affd_rest_controller::get_aff_structure_optimized_with_stats, "id", "asymlist");

		map_post_request("aff", &affd_rest_controller::post_custom_structure, "structure");
	}

	zh::reply get_aff_structure(const std::string &af_id);
	zeep::json::element get_aff_status(const std::string &af_id);

	zeep::json::element get_aff_structure_json(const std::string &af_id);
	zeep::json::element get_aff_3d_beacon(std::string id);

	zh::reply get_aff_structure_stripped_def(const std::string &id, const std::string &asyms)
	{
		return get_aff_structure_stripped(id, asyms, 0);
	}

	zh::reply get_aff_structure_stripped(const std::string &af_id, const std::string &asyms, int identity);

	zh::reply get_aff_structure_optimized(const std::string &af_id, const std::string &asyms);
	zh::reply get_aff_structure_optimized_with_stats(const std::string &af_id, const std::string &asyms);

	// --------------------------------------------------------------------

	zeep::json::element post_custom_structure(const std::string &data);
};

zeep::json::element affd_rest_controller::get_aff_status(const std::string &af_id)
{
	const auto &[type, id, chunkNr] = parse_af_id(af_id);

	zeep::json::element result;

	if (type == EntryType::Custom)
	{
		const auto &[status, progress] = data_service::instance().get_status(id);

		result["status"] = status;
		result["progress"] = progress;
	}
	else if (fs::exists(file_locator::get_structure_file(id, chunkNr)))
		result["status"] = CustomStatus::Finished;
	else
		result["status"] = CustomStatus::Unknown;

	return result;
}

zh::reply affd_rest_controller::get_aff_structure(const std::string &af_id)
{
	zeep::http::reply rep(zeep::http::ok, {1, 1});

	const auto &[type, id, chunkNr] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(type, id, chunkNr);

	if (not fs::exists(file))
		return zeep::http::reply(zeep::http::not_found, {1, 1});

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		rep.set_content(new std::ifstream(file, std::ios::binary), "text/plain");
		rep.set_header("content-encoding", "gzip");
	}
	else
	{
		std::ifstream is(file, std::ios::binary);
		if (not is.is_open())
			return zeep::http::reply(zeep::http::not_found, {1, 1});

		io::filtering_stream<io::input> in;
		in.push(io::gzip_decompressor());
		in.push(is);

		std::ostringstream os;
		io::copy(in, os);

		rep.set_content(os.str(), "text/plain");
	}

	rep.set_header("content-disposition", "attachement; filename = \"AF-" + id + "-F" + std::to_string(chunkNr) + "-model_v1.cif\"");

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_stripped(const std::string &af_id, const std::string &asyms, int identity)
{
	zeep::http::reply rep(zeep::http::ok, {1, 1});

	std::set<std::string> requestedAsyms;
	ba::split(requestedAsyms, asyms, ba::is_any_of(","));

	std::unique_ptr<std::iostream> s(new std::stringstream);
	io::filtering_stream<io::output> out;

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		out.push(io::gzip_compressor());
		rep.set_header("content-encoding", "gzip");
	}

	out.push(*s.get());

	stripCifFile(af_id, requestedAsyms, identity, out);

	rep.set_content(s.release(), "text/plain");
	// rep.set_header("content-disposition", "attachement; filename = \"AF-" + id + "-F1-model_v1.cif\"");

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_optimized(const std::string &af_id, const std::string &asyms)
{
	zeep::http::reply rep(zeep::http::ok, {1, 1});

	std::set<std::string> requestedAsyms;
	ba::split(requestedAsyms, asyms, ba::is_any_of(","));

	std::unique_ptr<std::iostream> s(new std::stringstream);
	io::filtering_stream<io::output> out;

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		out.push(io::gzip_compressor());
		rep.set_header("content-encoding", "gzip");
	}

	out.push(*s.get());

	optimizeWithYasara("/opt/yasara/yasara", af_id, requestedAsyms, out);

	rep.set_content(s.release(), "text/plain");

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_optimized_with_stats(const std::string &af_id, const std::string &asyms)
{
	zeep::http::reply rep(zeep::http::ok, {1, 1});

	std::set<std::string> requestedAsyms;
	ba::split(requestedAsyms, asyms, ba::is_any_of(","));

	std::ostringstream os;

	auto stats = optimizeWithYasara("/opt/yasara/yasara", af_id, requestedAsyms, os);

	stats["model"] = os.str();

	std::unique_ptr<std::iostream> s(new std::stringstream);
	io::filtering_stream<io::output> out;

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		out.push(io::gzip_compressor());
		rep.set_header("content-encoding", "gzip");
	}

	out.push(*s.get());

	out << stats;

	rep.set_content(s.release(), "application/json");

	return rep;
}

zeep::json::element affd_rest_controller::get_aff_structure_json(const std::string &af_id)
{
	const auto &[type, id, chunkNr] = parse_af_id(af_id);

	fs::path file = file_locator::get_metadata_file(type, id, chunkNr);

	if (not fs::exists(file))
		throw zeep::http::not_found;

	zeep::json::element result;

	std::ifstream is(file);
	parse_json(is, result);

	return result;
}

zeep::json::element affd_rest_controller::get_aff_3d_beacon(std::string af_id)
{
	if (ba::ends_with(af_id, ".json"))
		af_id.erase(af_id.begin() + af_id.length() - 5, af_id.end());

	const auto &[type, id, chunkNr] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(type, id, chunkNr);

	if (not fs::exists(file))
		throw zeep::http::not_found;

	zeep::json::element result{
		{"uniprot_entry",
			{{"ac", id}}}};

	using namespace std::chrono;

	auto ft = fs::last_write_time(file);
	auto sctp = time_point_cast<system_clock::duration>(ft - decltype(ft)::clock::now() + system_clock::now());
	std::time_t cft = system_clock::to_time_t(sctp);
	std::tm *tm = std::gmtime(&cft);

	std::stringstream ss;
	ss << std::put_time(tm, "%F");

	// get the chain length...

	cif::file cf(file);
	if (cf.empty())
		throw zeep::http::not_found;
	auto &db = cf.front();
	// auto &struct_ref = db["struct_ref"];
	auto &struct_ref_seq = db["struct_ref_seq"];

	int uniprot_start, uniprot_end;
	cif::tie(uniprot_start, uniprot_end) = struct_ref_seq.front().get("db_align_beg", "db_align_end");

	result["structures"].push_back({{"model_identifier", id},
		{"model_category", "DEEP-LEARNING"},
		{"model_url", "https://alphafill.eu/v1/aff/" + id},
		{"model_page_url", "https://alphafill.eu/model?id=" + id},
		{"model_format", "MMCIF"},
		{"provider", "AlphaFill"},
		{"created", ss.str()},
		{"sequence_identity", 1.0},
		{"coverage", 1.0},
		{"uniprot_start", uniprot_start},
		{"uniprot_end", uniprot_end}});

	return result;
}

// --------------------------------------------------------------------

zeep::json::element affd_rest_controller::post_custom_structure(const std::string &data)
{
	auto hash = zeep::encode_hex(zeep::sha1(data));

	auto &&[status, progress] = data_service::instance().get_status(hash);

	if (status == CustomStatus::Unknown)
	{
		data_service::instance().queue(data, hash);
		status = CustomStatus::Queued;
	}

	return {
		{ "id", "CS-" + hash },
		{ "status", status }
	};
}

// --------------------------------------------------------------------

int server_main(int argc, char *const argv[])
{
	using namespace std::literals;

	int result = 0;

	zeep::value_serializer<CustomStatus>::instance("CustomStatus")
		("unknown", CustomStatus::Unknown)
		("queued", CustomStatus::Queued)
		("running", CustomStatus::Running)
		("finished", CustomStatus::Finished)
		("error", CustomStatus::Error);

	auto &config = cfg::config::instance();

	fs::path dbDir = config.get<std::string>("db-dir");

	LigandsTable::init(config.get<std::string>("ligands"));

	std::vector<std::string> vConn;
	std::string db_user;
	for (std::string opt : {"db-host", "db-port", "db-dbname", "db-user", "db-password"})
	{
		if (not config.has(opt))
			continue;

		vConn.push_back(opt.substr(3) + "=" + config.get<std::string>(opt));
		if (opt == "db-user")
			db_user = config.get<std::string>(opt);
	}

	db_connection::init(ba::join(vConn, " "));

	// --------------------------------------------------------------------

	if (config.has("rebuild-db"))
		return data_service::rebuild(db_user, dbDir);

	if (config.operands().size() < 2)
	{
		std::cout << "No command specified, use of of start, stop, status or reload" << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------

	if (config.has("db-link-template"))
		s_af_object.set_template(config.get<std::string>("db-link-template"));

	std::string command = config.operands()[1];

	zh::daemon server([&]()
		{
		// auto sc = new zh::security_context(secret, ibs_user_service::instance());
		// sc->add_rule("/admin", { "ADMIN" });
		// sc->add_rule("/admin/**", { "ADMIN" });
		// sc->add_rule("/media,/media/**", { "ADMIN", "EDITOR" });
		// sc->add_rule("/", {});

		auto s = new zeep::http::server(/*sc*/);

		if (config.has("context"))
			s->set_context_name(config.get<std::string>("context"));

		s->add_error_handler(new db_error_handler());
		s->add_error_handler(new missing_entry_error_handler());

#ifndef NDEBUG
		s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#else
		s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#endif
		// s->add_controller(new zh::login_controller());
		// s->add_controller(new user_admin_rest_controller());

		s->add_controller(new affd_html_controller());
		s->add_controller(new affd_rest_controller());

		return s; },
		PACKAGE_NAME);

	if (command == "start")
	{
		std::string address = "127.0.0.1";
		if (config.has("address"))
			address = config.get<std::string>("address");

		unsigned short port = 10342;
		if (config.has("port"))
			port = config.get<unsigned short>("port");

		std::string user = "www-data";
		if (config.has("user"))
			user = config.get<std::string>("user");

		std::cout << "starting server at http://" << address << ':' << port << '/' << std::endl;

		if (config.has("no-daemon"))
			result = server.run_foreground(address, port);
		else
			result = server.start(address, port, 4, 4, user);
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

	return result;
}
