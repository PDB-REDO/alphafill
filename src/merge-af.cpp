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

#include <fstream>
#include <iomanip>

#include <cif++/CifUtils.hpp>
#include <cif++/Structure.hpp>
#include <cif++/Point.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/iostreams/filtering_stream.hpp>
// #include <boost/numeric/ublas/matrix.hpp>
#include <boost/program_options.hpp>

#include <zeep/json/element.hpp>
#include <zeep/json/parser.hpp>

#include "revision.hpp"
#include "utilities.hpp"
// #include "ligands.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

double alignRange(mmcif::Structure &a, mmcif::Structure &b, size_t seqIDFirst, size_t seqIDLast)
{
	using namespace cif::literals;

	std::vector<mmcif::Point> pa, pb;

	auto &db_a = a.datablock();
	auto &db_b = b.datablock();

	auto &atom_site_a = db_a["atom_site"];
	auto &atom_site_b = db_b["atom_site"];

	auto &ma_qa_a = db_a["ma_qa_metric_local"];
	auto &ma_qa_b = db_a["ma_qa_metric_local"];

	auto ra = atom_site_a.find("pdbx_sifts_xref_db_num"_key >= seqIDFirst and "pdbx_sifts_xref_db_num"_key <= seqIDLast and "label_atom_id"_key != "OXT");
	auto rb = atom_site_b.find("pdbx_sifts_xref_db_num"_key >= seqIDFirst and "pdbx_sifts_xref_db_num"_key <= seqIDLast and "label_atom_id"_key != "OXT");

	for (auto rai = ra.begin(), rbi = rb.begin(); rai != ra.end() and rbi != rb.end(); ++rai, ++rbi)
	{
		if (rai == ra.end() or rbi == rb.end())
			throw std::runtime_error("Unequal number of atoms in ranges");

		std::string atom_id_a, atom_id_b, comp_id_a, comp_id_b;
		float xa, ya, za, xb, yb, zb;
		int seq_id_a, seq_id_b;

		std::tie(atom_id_a, comp_id_a, seq_id_a, xa, ya, za) = rai->get<std::string,std::string,int,float,float,float>({"label_atom_id", "label_comp_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"});
		std::tie(atom_id_b, comp_id_b, seq_id_b, xb, yb, zb) = rbi->get<std::string,std::string,int,float,float,float>({"label_atom_id", "label_comp_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z"});

		if (atom_id_a != atom_id_b or comp_id_a != comp_id_b)
			throw std::runtime_error("Atoms are not equal in ranges");
		
		if (atom_id_a != "CA" and atom_id_a != "C" and atom_id_a != "N" and atom_id_a != "O")
			continue;
		
		if (ma_qa_a.find1<float>("label_seq_id"_key == seq_id_a, "metric_value") < 70 or
			ma_qa_b.find1<float>("label_seq_id"_key == seq_id_b, "metric_value") < 70)
			continue;

		pa.emplace_back(xa, ya, za);
		pb.emplace_back(xb, yb, zb);
	}

	if (pa.empty())
		throw std::runtime_error("No atoms to align");

	auto ta = CenterPoints(pa);

	if (cif::VERBOSE > 0)
		std::cerr << "translate A: " << -ta << std::endl;

	auto tb = CenterPoints(pb);

	if (cif::VERBOSE > 0)
		std::cerr << "translate B: " << -tb << std::endl;

	auto rotation = AlignPoints(pa, pb);

	if (cif::VERBOSE > 0)
	{
		const auto &[angle, axis] = mmcif::QuaternionToAngleAxis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	for (auto &pt : pb)
		pt.rotate(rotation);

	double result = RMSd(pa, pb);

	if (cif::VERBOSE > 0)
		std::cerr << "RMSd: " << result << std::endl;

	return result;
}

void Extend(mmcif::Structure &s, std::vector<mmcif::Structure> &chunks, size_t seqIDFirst, size_t seqIDLast)
{
	struct SInfo {
		SInfo(mmcif::Structure &s)
			: s(s) {}

		mmcif::Structure &s;

	};

	std::vector<SInfo> parts;

	for (auto &chunk : chunks)
	{
		auto struct_ref = chunk.datablock()["struct_ref"].front();

		if (struct_ref["pdbx_align_end"].as<size_t>() < seqIDFirst)
			continue;
		
		if (struct_ref["pdbx_align_begin"].as<size_t>() > seqIDLast)
			continue;
		
		parts.emplace_back(chunk);
	}

	std::cout << "Merging " << parts.size() << " chunks for " << seqIDFirst << "-" << seqIDLast << std::endl;

	for (size_t i = 0; i + 1 < parts.size(); ++i)
	{
		for (size_t j = i + 1; j < parts.size(); ++j)
		{
			double rmsd = alignRange(parts[i].s, parts[j].s, seqIDFirst, seqIDLast);

			std::cout << i << "/" << j << " -> " << rmsd << std::endl;
		}
	}

}

// --------------------------------------------------------------------

int a_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] uniprot-id [output-file]"s);

	visible_options.add_options()
		("af-id", po::value<std::string>(), "AlphaFold ID");

	po::options_description hidden_options("hidden options");

	po::positional_options_description p;
	p.add("af-id", 1);

	po::variables_map vm = load_options(argc, argv, visible_options, hidden_options, p, "alphafill.conf");

	// --------------------------------------------------------------------

	fs::path afDir = vm["af-dir"].as<std::string>();
	if (not fs::is_directory(afDir))
		throw std::runtime_error("AlphfaFill data directory does not exist");

	file_locator::init(afDir, vm["structure-name-pattern"].as<std::string>(), vm["metadata-name-pattern"].as<std::string>());

	// --------------------------------------------------------------------

	if (vm.count("af-id") == 0)
	{
		std::cout << "AlphaFold ID not specified" << std::endl;
		exit(1);
	}

	std::string af_id = vm["af-id"].as<std::string>();

	auto af_files = file_locator::get_all_structure_files(af_id);

	if (af_files.size() < 2)
		throw std::runtime_error("Expected more than one AlphaFold structure files");

	std::vector<mmcif::File> files(af_files.size());
	std::vector<mmcif::Structure> structures;

	for (size_t i = 0; i < af_files.size(); ++i)
	{
		files[i].load(af_files[i]);
		structures.emplace_back(files[i]);
	}

	// This should be detected of course, but for now...
	const int kChunkSize = 200;

	// validate anyway, just to be sure
	size_t seqIDFirst = 1;
	size_t seqIDLast;

	for (auto &s : structures)
	{
		auto struct_ref = s.datablock()["struct_ref"].front();
		if (struct_ref["pdbx_align_begin"].as<size_t>() != seqIDFirst)
			throw std::runtime_error("Unexpected pdbx_align_begin, is the chunk size correct?");
		seqIDFirst += kChunkSize;
		seqIDLast = struct_ref["pdbx_align_end"].as<size_t>();
	}

	mmcif::Structure result(files[0]);

	seqIDFirst = 1;
	while (seqIDFirst < seqIDLast)
	{
		size_t seqIDEnd = seqIDFirst + kChunkSize - 1;
		if (seqIDEnd > seqIDLast)
			seqIDEnd = seqIDLast;

		Extend(result, structures, seqIDFirst, seqIDEnd);

		seqIDFirst += kChunkSize;
	}

	return 0;
}
