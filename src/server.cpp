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

#include "server.hpp"
#include "alphafill.hpp"
#include "data-service.hpp"
#include "db-connection.hpp"
#include "ligands.hpp"
#include "main.hpp"
#include "structure.hpp"
#include "utilities.hpp"

#include <zeep/crypto.hpp>
#include <zeep/http/daemon.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/http/uri.hpp>
#include <zeep/json/parser.hpp>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;
namespace zh = zeep::http;

const std::array<uint32_t, 6> kIdentities{ 70, 60, 50, 40, 30, 25 };

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
					id.erase(4);

				auto url = m_template;
				std::string::size_type s;
				while ((s = url.find("${id}")) != std::string::npos)
					url.replace(s, 5, id);

				to_element(result, url);
			}
			catch (const std::exception &e)
			{
				std::cerr << "Error getting post by ID: " << e.what() << '\n';
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
		mount("browserconfig.xml", &affd_html_controller::handle_file);
		mount("alphafill.json.schema", &affd_html_controller::schema);
		mount("favicon.ico", &affd_html_controller::handle_file);

		mount("structure-table-page", &affd_html_controller::structures_table);

		mount("{manual,man,genindex}/", &affd_html_controller::handle_help_file);
		mount("_static/", &affd_html_controller::handle_file);
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
	void handle_help_file(const zh::request &request, const zh::scope &scope, zh::reply &reply);
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
	std::vector<uint8_t> pae;
	double paeMean;
	double paeSD;
	bool firstHit = false;
	bool firstTransplant = false;
	int hitCount = 1;
	int transplantCount = 1;

	template <typename Archive>
	void serialize(Archive &ar, unsigned long)
	{
		ar &zeep::make_nvp("compound_id", compound_id)             //
			& zeep::make_nvp("analogue_id", analogue_id)           //
			& zeep::make_nvp("pdb_id", pdb_id)                     //
			& zeep::make_nvp("identity", identity)                 //
			& zeep::make_nvp("global-rmsd", gRMSd)                 //
			& zeep::make_nvp("asym_id", asym_id)                   //
			& zeep::make_nvp("local-rmsd", lRMSd)                  //
			& zeep::make_nvp("pae", pae)                           //
			& zeep::make_nvp("pae-mean", paeMean)                  //
			& zeep::make_nvp("pae-sd", paeSD)                      //
			& zeep::make_nvp("clash-score", clashScore)            //
			& zeep::make_nvp("first-hit", firstHit)                //
			& zeep::make_nvp("first-transplant", firstTransplant)  //
			& zeep::make_nvp("hit-count", hitCount)                //
			& zeep::make_nvp("transplant-count", transplantCount); //
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

	const std::string af_id = request.get_parameter("id");
	auto status = data_service::instance().get_status(af_id);

	switch (status.status)
	{
		case CustomStatus::Unknown:
			if (data_service::instance().exists_in_afdb(af_id))
			{
				auto id = data_service::instance().queue_af_id(af_id);
				sub.put("hash", id);
				sub.put("status", status.status);
				get_template_processor().create_reply_from_template("wait", sub, reply);
				return;
			}

			throw missing_entry_error(af_id);

		case CustomStatus::Finished:
			break;

		case CustomStatus::Error:
		{
			const auto &[type, afId, chunkNr, version] = parse_af_id(af_id);
			std::ifstream errmsg(file_locator::get_error_file(type, afId, chunkNr, version));

			sub.put("error-id", af_id);
			sub.put("error", status.message.value_or("<< message is missing >>"));

			get_template_processor().create_reply_from_template("index", sub, reply);
			return;
		}

		default:
			sub.put("hash", af_id);
			sub.put("status", status.status);
			get_template_processor().create_reply_from_template("wait", sub, reply);
			return;
	}

	const auto &[type, afId, chunkNr, version] = parse_af_id(af_id);

	fs::path jsonFile = file_locator::get_metadata_file(type, afId, chunkNr, version);
	fs::path cifFile = file_locator::get_structure_file(type, afId, chunkNr, version);

	sub.put("af_id", af_id);
	sub.put("id", afId);
	sub.put("chunk", chunkNr);
	sub.put("type", type);
	sub.put("version", version);

	bool chunked = type == EntryType::AlphaFold and (chunkNr > 1 or fs::exists(file_locator::get_metadata_file(type, afId, 2, version)));

	sub.put("chunked", chunked);

	if (chunked)
	{
		auto allChunks = file_locator::get_all_structure_files(afId, version);
		json chunks;

		for (size_t i = 0; i < allChunks.size(); ++i)
			chunks.emplace_back(afId + "-F" + std::to_string(i + 1));

		sub.put("chunks", chunks);
	}

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
				if (hit["alignment"]["identity"].as<double>() * 100 < i)
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
		double hitIdentity = hit["alignment"]["identity"].as<double>();
		if (hitIdentity * 100 < identity)
			continue;

		++hit_nr;
		for (auto &transplant : hit["transplants"])
		{
			std::vector<uint8_t> pae;
			if (transplant["pae"].is_object() and transplant["pae"]["matrix"].is_array())
			{
				for (auto &row : transplant["pae"]["matrix"])
				{
					for (auto &f : row)
						pae.push_back(f.as<uint8_t>());
				}
			}

			transplants.emplace_back(transplant_info{
				transplant["compound_id"].as<std::string>(),
				transplant["analogue_id"].as<std::string>(),
				hit_nr,
				hit["pdb_id"].as<std::string>() + '.' + hit["pdb_asym_id"].as<std::string>(),
				hitIdentity,
				hit["global_rmsd"].as<double>(),
				transplant["asym_id"].as<std::string>(),
				transplant["clash"]["score"].as<double>(),
				transplant["local_rmsd"].as<double>(),
				std::move(pae),
				transplant["pae"]["mean"].as<double>(),
				transplant["pae"]["stddev"].as<double>() });
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
		sub.put("cif-db-name", file.front().name());
	}
	catch (const std::exception &e)
	{
		sub.put("title", e.what());
	}

	// TODO: These magic numbers should of course be configurable parameters
	// 11.43, 3.04 voor global, en 3.10 en 0.92 voor local
	sub.put("cutoff", json{
						  // {"global", {{"unreliable", 11.43}, {"suspect", 3.04}}},
						  { "local", { { "unreliable", 3.10 }, { "suspect", 0.92 } } },
						  { "tcs", { { "unreliable", 1.27 }, { "suspect", 0.64 } } } });

	get_template_processor().create_reply_from_template("model", sub, reply);
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

	std::string af_id = request.get_parameter("id");
	const auto &[type, afId, chunkNr, version] = parse_af_id(af_id);

	sub.put("af_id", af_id);
	sub.put("chunk", chunkNr);
	sub.put("asym_id", asymID);
	sub.put("version", version);

	try
	{
		using namespace cif::literals;

		auto cif = file_locator::get_structure_file(type, afId, chunkNr, version);
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

	bool chunked = chunkNr > 1 or fs::exists(file_locator::get_metadata_file(type, afId, 2, version));

	sub.put("chunked", chunked);

	if (chunked)
	{
		auto allChunks = file_locator::get_all_structure_files(afId, version);
		json chunks;

		for (size_t i = 0; i < allChunks.size(); ++i)
			chunks.emplace_back(afId + "-F" + std::to_string(i + 1));

		sub.put("chunks", chunks);
	}

	get_template_processor().create_reply_from_template("optimized", sub, reply);
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

void affd_html_controller::handle_help_file(const zh::request &request, const zh::scope &scope, zh::reply &reply)
{
	zh::scope sub(scope);

	fs::path file = scope["baseuri"].as<std::string>();
	file /= "index.html";

	sub.put("baseuri", file.string());

	return get_template_processor().create_reply_from_template(file, sub, reply);
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

		map_get_request("aff/3d-beacon/{id}", &affd_rest_controller::get_aff_3d_beacon, "id", "version");

		map_get_request("aff/{id}/stripped/{asymlist}", &affd_rest_controller::get_aff_structure_stripped_def, "id", "asymlist");
		map_get_request("aff/{id}/stripped/{asymlist}/{identity}", &affd_rest_controller::get_aff_structure_stripped, "id", "asymlist", "identity");

		map_get_request("aff/{id}/optimized/{asymlist}", &affd_rest_controller::get_aff_structure_optimized, "id", "asymlist");
		map_get_request("aff/{id}/optimized-with-stats/{asymlist}", &affd_rest_controller::get_aff_structure_optimized_with_stats, "id", "asymlist");

		map_post_request("aff", &affd_rest_controller::post_custom_structure, "structure", "pae");
	}

	zh::reply get_aff_structure(const std::string &af_id);
	status_reply get_aff_status(const std::string &af_id);

	zeep::json::element get_aff_structure_json(const std::string &af_id);
	zeep::json::element get_aff_3d_beacon(std::string id, std::string version);

	zh::reply get_aff_structure_stripped_def(const std::string &id, const std::optional<std::string> &asyms)
	{
		return get_aff_structure_stripped(id, asyms, 0);
	}

	zh::reply get_aff_structure_stripped(const std::string &af_id, const std::optional<std::string> &asyms, int identity);

	zh::reply get_aff_structure_optimized(const std::string &af_id, const std::string &asyms);
	zh::reply get_aff_structure_optimized_with_stats(const std::string &af_id, const std::string &asyms);

	// --------------------------------------------------------------------

	zeep::json::element post_custom_structure(const std::string &data, const std::optional<std::string> pae);
};

status_reply affd_rest_controller::get_aff_status(const std::string &af_id)
{
	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);

	status_reply result = data_service::instance().get_status(af_id);

	if (result.status == CustomStatus::Unknown and fs::exists(file_locator::get_structure_file(id, chunkNr, version)))
		result.status = CustomStatus::Finished;

	return result;
}

zh::reply affd_rest_controller::get_aff_structure(const std::string &af_id)
{
	zeep::http::reply rep(zeep::http::ok, { 1, 1 });

	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(type, id, chunkNr, version);

	if (not fs::exists(file))
		return zeep::http::reply(zeep::http::not_found, { 1, 1 });

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		rep.set_content(new std::ifstream(file, std::ios::binary), "text/plain");
		rep.set_header("content-encoding", "gzip");
	}
	else
	{
		cif::gzio::ifstream in(file);

		if (not in.is_open())
			return zeep::http::reply(zeep::http::not_found, { 1, 1 });

		std::ostringstream os;
		os << in.rdbuf();

		rep.set_content(os.str(), "text/plain");
	}

	std::string filename = af_id + ".cif";
	rep.set_header("content-disposition", "attachement; filename = \"" + filename + '\"');

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_stripped(const std::string &af_id, const std::optional<std::string> &asyms, int identity)
{
	zeep::http::reply rep(zeep::http::ok, { 1, 1 });

	std::set<std::string> requestedAsymSet;
	if (asyms.has_value())
	{
		auto requestedAsyms = cif::split(*asyms, ",");
		requestedAsymSet = { requestedAsyms.begin(), requestedAsyms.end() };
	}

	std::unique_ptr<std::iostream> s(new std::stringstream);

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		cif::gzio::basic_ogzip_streambuf<char, std::char_traits<char>> buffer;
		buffer.init(s->rdbuf());
		std::ostream os(&buffer);

		stripCifFile(af_id, requestedAsymSet, identity, os);

		buffer.close();

		rep.set_header("content-encoding", "gzip");
	}
	else
		stripCifFile(af_id, requestedAsymSet, identity, *s);

	rep.set_content(s.release(), "text/plain");

	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);
	fs::path file = file_locator::get_structure_file(type, id, chunkNr, version);
	if (fs::exists(file))
	{
		auto filename = file.filename();
		if (filename.extension() == ".gz")
			filename.replace_extension("");

		rep.set_header("content-disposition", "attachement; filename = \"" + filename.string() + "\"");
	}

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_optimized(const std::string &af_id, const std::string &asyms)
{
	zeep::http::reply rep(zeep::http::ok, { 1, 1 });

	auto requestedAsyms = cif::split(asyms, ",");
	std::set<std::string> requestedAsymSet{ requestedAsyms.begin(), requestedAsyms.end() };

	std::unique_ptr<std::iostream> s(new std::stringstream);

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		cif::gzio::basic_ogzip_streambuf<char, std::char_traits<char>> buffer;
		buffer.init(s->rdbuf());
		std::ostream os(&buffer);

		optimizeWithYasara(af_id, requestedAsymSet, os);

		rep.set_header("content-encoding", "gzip");
	}
	else
		optimizeWithYasara(af_id, requestedAsymSet, *s);

	rep.set_content(s.release(), "text/plain");

	return rep;
}

zh::reply affd_rest_controller::get_aff_structure_optimized_with_stats(const std::string &af_id, const std::string &asyms)
{
	zeep::http::reply rep(zeep::http::ok, { 1, 1 });

	auto requestedAsyms = cif::split(asyms, ",");
	std::set<std::string> requestedAsymSet{ requestedAsyms.begin(), requestedAsyms.end() };

	std::ostringstream os;

	auto stats = optimizeWithYasara(af_id, requestedAsymSet, os);

	stats["model"] = os.str();

	std::unique_ptr<std::iostream> s(new std::stringstream);

	if (get_header("accept-encoding").find("gzip") != std::string::npos)
	{
		cif::gzio::basic_ogzip_streambuf<char, std::char_traits<char>> buffer;
		buffer.init(s->rdbuf());
		std::ostream os(&buffer);

		os << stats;

		rep.set_header("content-encoding", "gzip");
	}
	else
		*s << stats;

	rep.set_content(s.release(), "application/json");

	return rep;
}

zeep::json::element affd_rest_controller::get_aff_structure_json(const std::string &af_id)
{
	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);

	fs::path file = file_locator::get_metadata_file(type, id, chunkNr, version);

	if (not fs::exists(file))
		throw zeep::http::not_found;

	zeep::json::element result;

	std::ifstream is(file);
	parse_json(is, result);

	return result;
}

zeep::json::element affd_rest_controller::get_aff_3d_beacon(std::string af_id, std::string version_3dbeacons)
{
	using namespace cif::literals;

	static const std::regex KVersionRX(R"((\d+)(?:\.(\d+))?(?:\.(\d+))?)");

	if (cif::ends_with(af_id, ".json"))
		af_id.erase(af_id.begin() + af_id.length() - 5, af_id.end());

	const auto &[type, id, chunkNr, version] = parse_af_id(af_id);

	fs::path file = file_locator::get_structure_file(type, id, chunkNr, version);

	if (not fs::exists(file))
	{
		// Disabled due to too many requests
		// auto &data_service = data_service::instance();
		// data_service.queue_3d_beacon_request(id);

		throw zeep::http::not_found;
	}

	int version_major = 1 /* , version_minor = 0 */;
	std::smatch m;
	if (std::regex_match(version_3dbeacons, m, KVersionRX))
	{
		version_major = std::stoi(m[1]);
		// if (m.length() >= 2)
		// 	version_minor = std::stoi(m[2]);
	}

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
	auto &struct_ref = db["struct_ref"];
	auto &struct_ref_seq = db["struct_ref_seq"];

	int uniprot_start, uniprot_end;
	cif::tie(uniprot_start, uniprot_end) = struct_ref_seq.front().get("db_align_beg", "db_align_end");

	std::string db_code = struct_ref.front()["db_code"].as<std::string>();

	zeep::json::element result{
		{ "uniprot_entry", { //
							   { "ac", id },
							   { "id", db_code },
							   { "sequence_length", uniprot_end - uniprot_start + 1 } } }
	};

	if (version_major >= 2)
	{
		zeep::json::element summary{
			{ "model_identifier", id },
			{ "model_category", "TEMPLATE-BASED" },
			{ "model_url", "https://alphafill.eu/v1/aff/" + id },
			{ "model_format", "MMCIF" },
			{ "model_page_url", "https://alphafill.eu/model?id=" + id },
			{ "provider", "AlphaFill" },
			{ "created", ss.str() },
			{ "sequence_identity", 1.0 },
			{ "uniprot_start", uniprot_start },
			{ "uniprot_end", uniprot_end },
			{ "coverage", 1.0 },
		};

		auto &entities = summary["entities"];
		auto &struct_asym = db["struct_asym"];

		for (const auto &[id, description, type] : db["entity"].rows<int, std::string, std::string>("id", "pdbx_description", "type"))
		{
			if (type == "polymer")
			{
				entities.push_back({ { "entity_type", "POLYMER" },
					{ "entity_poly_type", "POLYPEPTIDE(L)" },
					{ "description", description } });
				entities.back()["chain_ids"].push_back("A");
				continue;
			}

			if (type == "non-polymer")
			{
				entities.push_back({ { "entity_type", "NON-POLYMER" },
					{ "description", description } });

				auto &chain_ids = entities.back()["chain_ids"];

				for (auto asym_id : struct_asym.find<std::string>("entity_id"_key == id, "id"))
					chain_ids.push_back(asym_id);

				continue;
			}
		}

		result["structures"].push_back({ { "summary", summary } });
	}
	else
	{
		result["structures"].push_back({ { "model_identifier", id },
			{ "model_category", "DEEP-LEARNING" },
			{ "model_url", "https://alphafill.eu/v1/aff/" + id },
			{ "model_page_url", "https://alphafill.eu/model?id=" + id },
			{ "model_format", "MMCIF" },
			{ "provider", "AlphaFill" },
			{ "created", ss.str() },
			{ "sequence_identity", 1.0 },
			{ "coverage", 1.0 },
			{ "uniprot_start", uniprot_start },
			{ "uniprot_end", uniprot_end } });
	}

	return result;
}

// --------------------------------------------------------------------

zeep::json::element affd_rest_controller::post_custom_structure(const std::string &data, const std::optional<std::string> pae)
{
	auto id = "CS-" + zeep::encode_hex(zeep::sha1(data));

	auto status = data_service::instance().get_status(id);

	if (status.status == CustomStatus::Unknown)
	{
		data_service::instance().queue(data, pae, id);
		status.status = CustomStatus::Queued;
	}

	return {
		{ "id", id },
		{ "status", status.status }
	};
}

// --------------------------------------------------------------------

int server_main(int argc, char *const argv[])
{
	using namespace std::literals;

	auto &config = load_and_init_config(
		R"(usage: alphafill server <command> [options]

  where command is one of:

    start          Start a new webserver instance
    stop           Stop the currently running webserver instance
    status         Report the status of a currently running webserver
    reload         Reload the configuration, re-open log files and restart
                   the currently running webserver
)",

		mcfp::make_option<std::string>("db-dir", "Directory containing the alphafilled data"),
		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),

		mcfp::make_option<std::string>("pdb-fasta", "The FastA file containing the PDB sequences"),

		mcfp::make_option<std::string>("ligands", "af-ligands.cif", "File in CIF format describing the ligands and their modifications"),

		mcfp::make_option<float>("max-ligand-to-backbone-distance", 6, "The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure"),
		mcfp::make_option<float>("min-hsp-identity", 0.25, "The minimal identity for a high scoring pair (note, value between 0 and 1)"),
		mcfp::make_option<int>("min-alignment-length", 85, "The minimal length of an alignment"),
		mcfp::make_option<float>("min-separation-distance", 3.5, "The centroids of two identical ligands should be at least this far apart to count as separate occurrences"),
		mcfp::make_option<uint32_t>("blast-report-limit", 250, "Number of blast hits to use"),

		mcfp::make_option<std::string>("blast-matrix", "BLOSUM62", "Blast matrix to use"),
		mcfp::make_option<int>("blast-word-size", 3, "Blast word size"),
		mcfp::make_option<double>("blast-expect", 10, "Blast expect cut off"),
		mcfp::make_option("blast-no-filter", "Blast option for filter, default is to use low complexity filter"),
		mcfp::make_option("blast-no-gapped", "Blast option for gapped, default is to do gapped"),
		mcfp::make_option<int>("blast-gap-open", 11, "Blast penalty for gap open"),
		mcfp::make_option<int>("blast-gap-extend", 1, "Blast penalty for gap extend"),

		mcfp::make_option<float>("clash-distance-cutoff", 4, "The max distance between polymer atoms and ligand atoms used in calculating clash scores"),

		mcfp::make_option<float>("max-ligand-to-polymer-atom-distance", 6,
			"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure (for validation)"),

		mcfp::make_option<std::string>("structure-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.gz", "Pattern for locating structure files"),
		mcfp::make_option<std::string>("metadata-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.json", "Pattern for locating metadata files"),
		mcfp::make_option<std::string>("pdb-name-pattern", "${pdb-dir}/${id:1:2}/${id}/${id}_final.cif", "Pattern for locating PDB files"),

		mcfp::make_option<size_t>("threads,t", std::thread::hardware_concurrency(), "Number of threads to use, zero means all available cores"),

		mcfp::make_option<std::string>("alphafold-3d-beacon", "https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/${id}.json?provider=alphafold",
			"The URL of the 3d-beacons service for alphafold"),

		mcfp::make_option("no-daemon,F", "Do not fork a background process"),
		mcfp::make_option<std::string>("address", "127.0.0.1", "Address to listen to"),
		mcfp::make_option<unsigned short>("port", 10342, "Port to listen to"),
		mcfp::make_option<std::string>("user", "www-data", "User to run as"),
		mcfp::make_option<std::string>("context", "Reverse proxy context"),
		mcfp::make_option<std::string>("db-link-template", "Template for links to pdb(-redo) entry"),

		mcfp::make_option<std::string>("db-dbname", "AF DB name"),
		mcfp::make_option<std::string>("db-user", "AF DB owner"),
		mcfp::make_option<std::string>("db-password", "AF DB password"),
		mcfp::make_option<std::string>("db-host", "AF DB host"),
		mcfp::make_option<std::string>("db-port", "AF DB port"),

		mcfp::make_option<std::string>("custom-dir", (fs::temp_directory_path() / "alphafill").string(), "Directory for custom built entries"),

		mcfp::make_option<std::string>("yasara", "/opt/yasara/yasara", "Location of the yasara executable, needed for optimising")

	);

	parse_argv(argc, argv, config);

	// --------------------------------------------------------------------

	check_blast_index();

	if (not config.has("db-dir"))
	{
		std::cout << "Data directory not specified\n";
		return 1;
	}

	int result = 0;

	zeep::value_serializer<CustomStatus>::instance("CustomStatus") //
		("unknown", CustomStatus::Unknown)                         //
		("queued", CustomStatus::Queued)                           //
		("running", CustomStatus::Running)                         //
		("finished", CustomStatus::Finished)                       //
		("error", CustomStatus::Error);

	zeep::value_serializer<EntryType>::instance("EntryType") //
		("unknown", EntryType::Unknown)                      //
		("alphafold", EntryType::AlphaFold)                  //
		("custom", EntryType::Custom);

	LigandsTable::init(config.get("ligands"));

	std::vector<std::string> vConn;
	std::string db_user;
	for (std::string opt : { "db-host", "db-port", "db-dbname", "db-user", "db-password" })
	{
		if (not config.has(opt))
			continue;

		vConn.push_back(opt.substr(3) + "=" + config.get(opt));
		if (opt == "db-user")
			db_user = config.get(opt);
	}

	db_connection::init(cif::join(vConn, " "));

	// --------------------------------------------------------------------

	if (config.operands().size() != 1)
	{
		std::cout << "No command specified, use of of start, stop, status or reload\n";
		exit(1);
	}

	// --------------------------------------------------------------------

	if (config.has("db-link-template"))
		s_af_object.set_template(config.get("db-link-template"));

	std::string command = config.operands().front();

	zh::daemon server([&, nr_of_threads = config.get<size_t>("threads")]() mutable
		{

		auto &ds = data_service::instance();

		if (nr_of_threads == 0)
			nr_of_threads = std::thread::hardware_concurrency();

		ds.start_queue(nr_of_threads);

		auto s = new zeep::http::server(/*sc*/);

		if (config.has("context"))
			s->set_context_name(config.get("context"));

		s->add_error_handler(new db_error_handler());
		s->add_error_handler(new missing_entry_error_handler());

#if not defined(NDEBUG)
		s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#elif defined(WEBAPP_USES_RESOURCES) and WEBAPP_USES_RESOURCES
		s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#else
		s->set_template_processor(new zeep::http::file_based_html_template_processor(ALPHAFILL_DATA_DIR "/docroot"));
#endif

		s->add_controller(new affd_html_controller());
		s->add_controller(new affd_rest_controller());

		return s; },
		"alphafill");

	if (command == "start")
	{
		std::string address = config.get("address");
		unsigned short port = config.get<unsigned short>("port");
		std::string user = config.get("user");

		std::cout << "starting server at http://" << address << ':' << port << '/' << '\n';

		if (config.has("no-daemon"))
			result = server.run_foreground(address, port);
		else
			result = server.start(address, port, 1, 16, user);
	}
	else if (command == "stop")
		result = server.stop();
	else if (command == "status")
		result = server.status();
	else if (command == "reload")
		result = server.reload();
	else
	{
		std::cerr << "Invalid command " << std::quoted(command) << "\n";
		result = 1;
	}

	return result;
}
