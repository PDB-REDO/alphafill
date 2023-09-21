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

#include <iostream>
#include <filesystem>
#include <thread>

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>

#include "alphafill.hpp"
#include "data-service.hpp"
#include "main.hpp"
#include "revision.hpp"
#include "server.hpp"
#include "validate.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

int test_main(int argc, char *const argv[])
{
	data_service::instance().get_pae("P29376", 1, 3);

	return 0;
}


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

void parse_argv(int argc, char * const argv[], mcfp::config &config)
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

	config.parse_config_file("config", "alphafill.conf", { fs::current_path().string(), "/etc/" });
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

		auto &config = load_and_init_config(
R"(usage: alphafill command [options]

  where command is one of

    create-index   Create a FastA file based on data in the PDB files
                   (A FastA file is required to process files)
    process        Process an AlphaFill structure
    rebuild-db     Rebuild the databank
    server         Start a web server instance
)");

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
				std::cerr << "Missing command" << "\n\n";

			std::cerr << config << '\n';
			return config.has("help") ? 0 : 1;
		}

		// --------------------------------------------------------------------

		std::string command = config.operands().front();

		if (command == "create-index")
			result = create_index(argc - 1, argv + 1);
		else if (command == "process")
			result = alphafill_main(argc - 1, argv + 1);
		else if (command == "rebuild-db")
			result = rebuild_db_main(argc - 1, argv + 1);
		else if (command == "server")
			result = server_main(argc - 1, argv + 1);
		else 
		{
			std::cerr << "Unknown command\n\n"
					  << config << '\n';
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
