/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
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

#include "alphafill.hpp"
#include "config.hpp"
#include "main.hpp"
#include "revision.hpp"
#include "validate.hpp"

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include <filesystem>
#include <iostream>
#include <thread>

#if defined(BUILD_WEB_APPLICATION)
#include "data-service.hpp"
#include "db-connection.hpp"
#include "server.hpp"
#endif

namespace fs = std::filesystem;

// --------------------------------------------------------------------

int test_main(int argc, char *const argv[])
{
	return 0;
}

// --------------------------------------------------------------------

#if defined(BUILD_WEB_APPLICATION)
int rebuild_db_main(int argc, char *const argv[])
{
	using namespace std::literals;

	auto &config = load_and_init_config(
		"usage: alphafill rebuild-db [options]",
		mcfp::make_option<std::string>("db-dir", "Directory containing the alphafilled data"),
		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),

		mcfp::make_option<std::string>("structure-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.gz", "Pattern for locating structure files"),
		mcfp::make_option<std::string>("metadata-name-pattern", "${db-dir}/${id:0:2}/AF-${id}-F${chunk}-filled_v${version}.cif.json", "Pattern for locating metadata files"),
		mcfp::make_option<std::string>("pdb-name-pattern", "${pdb-dir}/${id:1:2}/${id}/${id}_final.cif", "Pattern for locating PDB files"),

		mcfp::make_option<std::string>("db-dbname", "AF DB name"),
		mcfp::make_option<std::string>("db-user", "AF DB owner"),
		mcfp::make_option<std::string>("db-password", "AF DB password"),
		mcfp::make_option<std::string>("db-host", "AF DB host"),
		mcfp::make_option<std::string>("db-port", "AF DB port"),

		mcfp::make_option<size_t>("threads,t", std::thread::hardware_concurrency(), "Number of threads to use, zero means all available cores"),

		mcfp::make_hidden_option<std::string>("custom-dir", (fs::temp_directory_path() / "alphafill").string(), "Directory for custom built entries")
		);

	parse_argv(argc, argv, config);

	// --------------------------------------------------------------------

	if (not config.has("db-dir"))
	{
		std::cout << "Data directory not specified\n";
		return 1;
	}

	// --------------------------------------------------------------------

	fs::path dbDir = config.get("db-dir");

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

	return data_service::rebuild(db_user, dbDir);
}
#endif

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << '\n';
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

void parse_argv(int argc, char *const argv[], mcfp::config &config)
{
	config.parse(argc, argv);

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help"))
	{
		std::cerr << config << '\n';
		exit(config.has("help") ? 0 : 1);
	}

	if (config.has("quiet"))
		cif::VERBOSE = -1;
	else
		cif::VERBOSE = config.count("verbose");

	config.set_ignore_unknown(true);

	std::error_code ec;
	config.parse_config_file("config", "alphafill.conf", { fs::current_path().string(), CONFIG_PATH }, ec);
	if (ec and ec != mcfp::config_error::config_file_not_found)
		throw std::system_error(ec, "Error parsing config file");
}

// --------------------------------------------------------------------

int main(int argc, char *const argv[])
{
	int result = 0;

	try
	{
#if defined(ALPHAFILL_DATA_DIR)
		cif::add_data_directory(ALPHAFILL_DATA_DIR);
#endif

		std::string command;
		const char *exe = strrchr(argv[0], '/');
		if (exe == nullptr)
			exe = argv[0];
		else
			++exe;

		if (exe != nullptr and strncmp(exe, "alphafill-", 10) == 0)
		{
			auto &config = load_and_init_config(argv[0]);

			std::error_code ec;
			config.set_ignore_unknown(true);
			config.parse(argc, argv, ec);

			if (config.has("version"))
			{
				write_version_string(std::cout, config.has("verbose"));
				exit(0);
			}

			command = exe + 10;
		}
		else
		{
			const std::string usage =
R"(usage: alphafill command [options]

where command is one of

    create-index   Create a FastA file based on data in the PDB files
                   (A FastA file is required to process files)
    process        Process an AlphaFill structure)"
#if defined(BUILD_WEB_APPLICATION)
R"(
    rebuild-db     Rebuild the databank
    server         Start a web server instance)"
#endif
R"(

The following options are always recognized:
)";

			auto &config = load_and_init_config(usage);

			std::error_code ec;
			config.set_ignore_unknown(true);
			config.parse(argc, argv, ec);

			if (config.has("version"))
			{
				write_version_string(std::cout, config.has("verbose"));
				exit(0);
			}

			if (config.operands().empty() or ec)
			{
				if (ec)
					std::cerr << "Error parsing arguments: " << ec.message() << "\n\n";

				if (config.operands().empty())
					std::cerr << "Missing command"
							<< "\n\n";

				std::cerr << config << '\n';
				return config.has("help") ? 0 : 1;
			}

			command = config.operands().front();
		}

		// --------------------------------------------------------------------

		if (command == "create-index")
			result = create_index(argc - 1, argv + 1);
		else if (command == "process")
			result = alphafill_main(argc - 1, argv + 1);
#if defined(BUILD_WEB_APPLICATION)
		else if (command == "rebuild-db")
			result = rebuild_db_main(argc - 1, argv + 1);
		else if (command == "server")
			result = server_main(argc - 1, argv + 1);
#endif
		else
		{
			std::cerr << "Unknown command " << std::quoted(command) << "\n\n"
					  << mcfp::config::instance() << '\n';
			result = 1;
		}
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		result = 1;
	}

	return result;
}
