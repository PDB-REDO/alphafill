/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>

#include <zeep/json/element.hpp>

#include "blast.hpp"
#include "ligands.hpp"
#include "queue.hpp"
#include "revision.hpp"
#include "utilities.hpp"
#include "validate.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

using json = zeep::json::element;

// --------------------------------------------------------------------

fs::path pdbFileForID(const fs::path &pdbDir, std::string pdb_id)
{
	for (auto &ch : pdb_id)
		ch = std::tolower(ch);

	// try a PDB-REDO layout first
	fs::path pdb_path = pdbDir / pdb_id.substr(1, 2) / pdb_id / (pdb_id + "_final.cif");
	if (not fs::exists(pdb_path))
		pdb_path = pdbDir / pdb_id.substr(1, 2) / (pdb_id + ".cif.gz");

	if (not fs::exists(pdb_path))
		throw std::runtime_error("PDB file for " + pdb_id + " not found");

	return pdb_path;
}

sequence getSequenceForStrand(cif::Datablock &db, const std::string &strand)
{
	using namespace cif::literals;

	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];
	auto r = pdbx_poly_seq_scheme.find("pdb_strand_id"_key == strand);
	if (r.empty())
		throw std::runtime_error("Could not locate sequence in PDB for strand id " + strand);

	auto entity_id = r.front()["entity_id"].as<std::string>();

	auto &entity_poly = db["entity_poly"];
	auto pdb_seq = entity_poly.find1<std::string>("entity_id"_key == entity_id, "pdbx_seq_one_letter_code_can");

	return encode(pdb_seq);
}

int validateFastA(fs::path fasta, fs::path dbDir, int threads)
{
	int result = 0;

	std::ifstream f(fasta);
	if (not f.is_open())
		throw std::runtime_error("Could not open fasta file");
	
	std::string id, strand, seq, line;

	std::vector<std::string> mismatch, unequal_length, not_x_related;
	blocking_queue<std::tuple<std::string,std::string,std::string>> q;

	std::mutex guard;

	std::vector<std::thread> t;
	for (int i = 0; i < threads; ++i)
	{
		t.emplace_back([&]()
		{
			for (;;)
			{
				const auto &&[id, strand, seq] = q.pop();

				if (id.empty())	// sentinel
				{
					q.push({});
					break;
				}

				try
				{
					cif::File pdbFile(pdbFileForID(dbDir, id));

					auto a = getSequenceForStrand(pdbFile.firstDatablock(), strand);
					auto b = encode(seq);

					if (a == b)
						continue;

					std::unique_lock lock(guard);

					result = -1;

					std::cerr << "Mismatch for " << id << " strand " << strand << std::endl;

					std::cerr << std::endl
								<< decode(a) << std::endl
								<< seq << std::endl
								<< std::endl;

					mismatch.push_back(id);

					if (a.length() != b.length())
					{
						unequal_length.push_back(id);
						continue;
					}

					for (size_t i = 0; i < a.length(); ++i)
					{
						if (a[i] == b[i])
							continue;
						
						if (a[i] == 22 or b[i] == 22)	// either one is X
							continue;
						
						not_x_related.push_back(id);
					}
				}
				catch (std::exception const &ex)
				{
					std::cerr << ex.what() << std::endl;
				}
			}
		});
	}

	while (std::getline(f, line))
	{
		if (line[0] == '>')
		{
			if (not id.empty())
			{
				if (seq.empty() or strand.empty())
					throw std::runtime_error("No sequence for id " + id);

				q.push({id, strand, seq });
			}
				
			id = line.substr(1, 4);
			strand = line.substr(6);
			seq.clear();
		}
		else
			seq.append(line.begin(), line.end());
	}

	// signal end
	q.push({});

	for (auto &ti: t)
		ti.join();

	if (result)
	{
		mismatch.erase(std::unique(mismatch.begin(), mismatch.end()), mismatch.end());
		unequal_length.erase(std::unique(unequal_length.begin(), unequal_length.end()), unequal_length.end());
		not_x_related.erase(std::unique(not_x_related.begin(), not_x_related.end()), not_x_related.end());

		std::cout << "Report for fasta check" << std::endl
				  << std::string(80, '-') << std::endl
				  << std::endl
				  << "PDB ID's with mismatches" << std::endl
				  << ba::join(mismatch, ", ") << std::endl
				  << std::endl
				  << "PDB ID's with differing sequence length" << std::endl
				  << ba::join(unequal_length, ", ") << std::endl
				  << std::endl
				  << "PDB ID's with mismatches that do not involve X" << std::endl
				  << ba::join(not_x_related, ", ") << std::endl
				  << std::endl;
	}

	return result;
}

// --------------------------------------------------------------------

using mmcif::Point;
using mmcif::Quaternion;

bool validateHit(mmcif::Structure &structure, const BlastHit &hit)
{
	using namespace cif::literals;

	return getSequenceForStrand(structure.datablock(), hit.mDefLine.substr(6)) == hit.mTarget;
}

std::vector<mmcif::Residue *> getResiduesForChain(mmcif::Structure &structure, const std::string &chain_id)
{
	std::vector<mmcif::Residue *> result;

	for (auto &poly : structure.polymers())
	{
		if (poly.chainID() != chain_id)
			continue;

		for (auto &res : poly)
			result.emplace_back(&res);
	}

	return result;
}

std::tuple<std::vector<size_t>, std::vector<size_t>> getTrimmedIndicesForHsp(const BlastHsp &hsp)
{
	std::vector<size_t> ixq, ixt;

	assert(hsp.mAlignedQuery.length() == hsp.mAlignedTarget.length());

	auto &qa = hsp.mAlignedQuery;
	auto &ta = hsp.mAlignedTarget;

	assert(qa.length() == ta.length());

	size_t qix = hsp.mQueryStart, tix = hsp.mTargetStart;

	for (size_t i = 0; i < qa.length(); ++i)
	{
		if (is_gap(qa[i]))
		{
			++tix;
			continue;
		}

		if (is_gap(ta[i]))
		{
			++qix;
			continue;
		}

		ixq.push_back(qix++);
		ixt.push_back(tix++);
	}

	assert(qix == hsp.mQueryEnd);
	assert(tix == hsp.mTargetEnd);

	return {ixq, ixt};
}

std::tuple<std::vector<Point>, std::vector<Point>> selectAtomsNearResidue(
	const std::vector<mmcif::Residue *> &pdb, const std::vector<size_t> &pdb_ix,
	const std::vector<mmcif::Residue *> &af, const std::vector<size_t> &af_ix,
	const std::vector<mmcif::Atom> &residue, float maxDistance, const Ligand &ligand)
{
	std::vector<Point> ra, rb;

	float maxDistanceSq = maxDistance * maxDistance;

	assert(pdb_ix.size() == af_ix.size());

	for (size_t i = 0; i < pdb_ix.size(); ++i)
	{
		bool nearby = false;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(pdb_ix[i] < pdb.size());

			auto atom = pdb[pdb_ix[i]]->atomByID(atom_id);
			if (not atom)
				continue;

			for (auto &b : residue)
			{
				if (ligand.drops(b.labelAtomID()))
					continue;

				if (DistanceSquared(atom, b) > maxDistanceSq)
					continue;

				nearby = true;
				break;
			}

			if (nearby)
				break;
		}

		if (not nearby)
			continue;

		for (const char *atom_id : {"C", "CA", "N", "O"})
		{
			assert(af_ix[i] < af.size());
			assert(pdb_ix[i] < pdb.size());

			auto pt_a = pdb[pdb_ix[i]]->atomByID(atom_id);
			auto pt_b = af[af_ix[i]]->atomByID(atom_id);

			if (not pt_a and pt_b)
				continue;

			ra.push_back(pt_a.location());
			rb.push_back(pt_b.location());
		}
	}

	return {ra, rb};
}

// --------------------------------------------------------------------

enum class UniqueType
{
	Seen, Unique, MoreAtoms
};

std::tuple<UniqueType,std::string> isUniqueLigand(const mmcif::Structure &structure, float minDistance, const mmcif::Residue &lig, std::string_view id)
{
	std::tuple<UniqueType,std::string> result{ UniqueType::Unique, "" };

	auto minDistanceSq = minDistance * minDistance;

	std::vector<mmcif::Point> pa;

	for (auto &a : lig.atoms())
		pa.push_back(a.location());
	auto ca = mmcif::Centroid(pa);

	for (auto &np : structure.nonPolymers())
	{
		if (np.compoundID() != id)
			continue;
		std::vector<mmcif::Point> pb;

		for (auto &a : np.atoms())
			pb.push_back(a.location());
		auto cb = mmcif::Centroid(pb);

		if (DistanceSquared(ca, cb) < minDistanceSq)
		{
			if (lig.unique_atoms().size() > np.unique_atoms().size())
				result = { UniqueType::MoreAtoms, np.asymID() };
			else
				result = { UniqueType::Seen, np.asymID() };

			break;
		}
	}

	return result;
}

// --------------------------------------------------------------------

int GeneratePDBList(fs::path pdbDir, LigandsTable &ligands, const std::string &output)
{
	std::cerr << "collecting files ";

	std::vector<fs::path> files;

	for (fs::directory_iterator iter(pdbDir); iter != fs::directory_iterator(); ++iter)
	{
		if (not iter->is_directory())
			continue;
		
		if (iter->path().filename().string().length() != 2)
			continue;
		
		// Regular PDB layout
		for (fs::directory_iterator fiter(iter->path()); fiter != fs::directory_iterator(); ++fiter)
		{
			if (not fiter->is_regular_file())
				continue;

			fs::path file = fiter->path();

			std::string name = file.filename().string();
			if (not ba::ends_with(name, ".cif.gz"))
				continue;

			files.push_back(file);

			if (files.size() % 10000 == 0)
				std::cerr << '.';
		}

		// PDB-REDO layout
		for (fs::recursive_directory_iterator fiter(iter->path()); fiter != fs::recursive_directory_iterator(); ++fiter)
		{
			fs::path file = fiter->path();

			std::string name = file.filename().string();
			if (not ba::ends_with(name, "_final.cif"))
				continue;

			files.push_back(file);

			if (files.size() % 10000 == 0)
				std::cerr << '.';
		}
	}

	std::cerr << " done!" << std::endl
			  << "need to process " << files.size() << " files" << std::endl;

	std::sort(files.begin(), files.end());

	cif::Progress progress(files.size(), "Parsing PDB files");

	std::vector<std::string> result;

	int nrOfThreads = std::thread::hardware_concurrency();
	blocking_queue<fs::path> q;
	std::mutex guard;

	std::vector<std::thread> t;
	for (int i = 0; i < nrOfThreads; ++i)
	{
		t.emplace_back([&]()
		{
			for (;;)
			{
				fs::path f = q.pop();

				if (f.empty())	// sentinel
				{
					q.push({});
					break;
				}

				try
				{
					cif::File file(f);

					progress.consumed(1);

					auto &db = file.firstDatablock();
					auto pdb_chem_comp = db.get("chem_comp");
					if (not pdb_chem_comp)
						continue;

					bool none = true;
					for (const auto &[comp_id] : pdb_chem_comp->rows<std::string>("id"))
					{
						if (ligands[comp_id])
						{
							none = false;
							break;
						}
					}

					if (none)
						continue;

					std::unique_lock lock(guard);
					result.push_back(db.getName());
				}
				catch(const std::exception& e)
				{
					std::cerr << e.what() << std::endl;
				}
			}
		});
	}

	for (auto &file : files)
		q.push(file);

	q.push({});

	// signal end

	for (auto &ti: t)
		ti.join();

	if (not output.empty())
	{
		std::ofstream outfile(output);
		for (auto &id : result)
			outfile << id << std::endl;
	}
	else
		std::cout << ba::join(result, ", ") << std::endl;

	return 0;
}

// --------------------------------------------------------------------

int a_main(int argc, char *const argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);

	visible_options.add_options()
		("pdb-fasta", po::value<std::string>(), "The FastA file containing the PDB sequences")
		("pdb-id-list", po::value<std::string>(), "Optional file containing the list of PDB ID's that have any of the transplantable ligands")

		("max-ligand-to-backbone-distance", po::value<float>()->default_value(6), "The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure")
		("min-hsp-identity", po::value<float>()->default_value(0.25), "The minimal identity for a high scoring pair (note, value between 0 and 1)")
		("min-alignment-length", po::value<int>()->default_value(85), "The minimal length of an alignment")
		("min-separation-distance", po::value<float>()->default_value(3.5), "The centroids of two identical ligands should be at least this far apart to count as separate occurrences")
		("blast-report-limit", po::value<uint32_t>()->default_value(250), "Number of blast hits to use")

		("clash-distance-cutoff", po::value<float>()->default_value(4), "The max distance between polymer atoms and ligand atoms used in calculating clash scores")

		("threads,t", po::value<size_t>()->default_value(1), "Number of threads to use, zero means all available cores");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("xyzin,i", po::value<std::string>(), "coordinates file")
		("output,o", po::value<std::string>(), "Output to this file")

		("validate-fasta", "Validate the FastA file (check if all sequence therein are the same as in the corresponding PDB files)")
		("prepare-pdb-list", "Generate a list with PDB ID's that contain any of the ligands")

		("test-pdb-id", po::value<std::string>(), "Test with single PDB ID")

		("test", "Run test code");

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);

	po::variables_map vm = load_options(argc, argv, visible_options, hidden_options, p, "alphafill.conf");

	// --------------------------------------------------------------------

	if (vm.count("pdb-fasta") == 0)
	{
		std::cout << "fasta file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("pdb-dir") == 0)
	{
		std::cout << "PDB directory not specified" << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------

	std::string fasta = vm["pdb-fasta"].as<std::string>();
	if (not fs::exists(fasta))
		throw std::runtime_error("PDB-Fasta file does not exist (" + fasta + ")");

	fs::path pdbDir = vm["pdb-dir"].as<std::string>();
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	// --------------------------------------------------------------------
	
	if (vm.count("validate-fasta"))
		return validateFastA(vm["pdb-fasta"].as<std::string>(), vm["pdb-dir"].as<std::string>(), std::thread::hardware_concurrency());

	// --------------------------------------------------------------------

	fs::path ligandsFile = vm["ligands"].as<std::string>();
	if (not fs::exists(ligandsFile))
	{
		std::cerr << "Ligands file not found" << std::endl;
		exit(1);
	}

	LigandsTable ligands(ligandsFile);

	float maxDistance = vm["clash-distance-cutoff"].as<float>();

	// --------------------------------------------------------------------
	
	if (vm.count("prepare-pdb-list"))
		return GeneratePDBList(vm["pdb-dir"].as<std::string>(), ligands, vm.count("output") ? vm["output"].as<std::string>() : "");

	cif::iset pdbIDsContainingLigands;
	if (vm.count("test-pdb-id"))
		pdbIDsContainingLigands.emplace(vm["test-pdb-id"].as<std::string>());
	else if (vm.count("pdb-id-list"))
	{
		std::ifstream file(vm["pdb-id-list"].as<std::string>());

		if (not file.is_open())
		{
			std::cerr << "Could not open pdb-id-list " << vm["pdb-id-list"].as<std::string>() << std::endl;
			exit(1);
		}

		std::string line;
		while (std::getline(file, line))
			pdbIDsContainingLigands.insert(line);
	}

	// --------------------------------------------------------------------
	
	if (vm.count("xyzin") == 0)
	{
		std::cout << "Input file not specified" << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------

	float minHspIdentity = vm["min-hsp-identity"].as<float>();
	size_t minAlignmentLength = vm["min-alignment-length"].as<int>();
	float minSeparationDistance = vm["min-separation-distance"].as<float>();
	float maxLigandBackboneDistance = vm["max-ligand-to-backbone-distance"].as<float>();

	// --------------------------------------------------------------------

	fs::path xyzin = vm["xyzin"].as<std::string>();
	mmcif::File f(xyzin);
	mmcif::Structure af_structure(f, 1, mmcif::StructureOpenOptions::SkipHydrogen);

	// --------------------------------------------------------------------
	// atoms for the clash score calculation
	std::vector<CAtom> polyAtoms;

	for (auto &res : af_structure.polymers().front())
	{
		for (auto &atom : res.atoms())
			polyAtoms.emplace_back(atom);
	}

	// --------------------------------------------------------------------
	// fetch the (single) chain

	const std::regex kIDRx(R"(^>(\w{4,8})_(\w)( .*)?)");

	auto &db = f.data();

	std::string afID;
	cif::tie(afID) = db["entry"].front().get("id");

	auto now = boost::posix_time::second_clock::universal_time();

	json result = {
		{"id", afID},
		{"file", xyzin.string()},
		{"date", to_iso_extended_string(now.date())},
		{"alphafill_version", kVersionNumber}};

	json &hits = result["hits"] = json::array();

	// keep a LRU cache of mmCIF parsed files
	std::list<std::tuple<std::string,std::shared_ptr<mmcif::File>>> mmCifFiles;

	for (auto r : db["entity_poly"])
	{
		auto &&[id, seq] = r.get<std::string, std::string>({"entity_id", "pdbx_seq_one_letter_code"});

		if (cif::VERBOSE > 0)
			std::cerr << "Blasting:" << std::endl
					  << seq << std::endl
					  << std::endl;

		size_t threads = vm["threads"].as<size_t>();
		if (threads == 0)
			threads = std::thread::hardware_concurrency();

		auto result = BlastP(fasta, seq, "BLOSUM62", 3, 10, true, true, 11, 1, vm["blast-report-limit"].as<uint32_t>(), threads);

		if (cif::VERBOSE > 0)
			std::cerr << "Found " << result.size() << " hits" << std::endl;

		std::unique_ptr<cif::Progress> progress;
		if (cif::VERBOSE < 1)
			progress.reset(new cif::Progress(result.size() + 1, "matching"));

		for (auto &hit : result)
		{
			if (progress)
				progress->consumed(1);

			std::smatch m;
			if (not regex_match(hit.mDefLine, m, kIDRx))
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Could not interpret defline from blast hit:" << std::endl
							  << hit.mDefLine << std::endl;
				continue;
			}

			std::string pdb_id = m[1].str();
			std::string chain_id = m[2].str();

			if (progress)
				progress->message(pdb_id);

			if (not (pdbIDsContainingLigands.empty() or pdbIDsContainingLigands.count(pdb_id)))
				continue;

			if (cif::VERBOSE > 0)
				std::cerr << "pdb id: " << pdb_id << '\t' << "chain id: " << chain_id << std::endl;
			
			try
			{
				auto ci = mmCifFiles.begin();

				for (; ci != mmCifFiles.end(); ++ci)
				{
					auto &&[c_id, c_file] = *ci;

					if (c_id != pdb_id)
						continue;
					
					if (ci == mmCifFiles.begin())
						break;
					
					mmCifFiles.emplace_front(c_id, c_file);
					mmCifFiles.erase(ci);

					ci = mmCifFiles.begin();
					break;
				}

				if (ci == mmCifFiles.end())
				{
					fs::path pdb_path = pdbFileForID(pdbDir, pdb_id);

					mmCifFiles.emplace_front(pdb_id, std::make_shared<mmcif::File>(pdb_path.string()));
					ci = mmCifFiles.begin();

					if (mmCifFiles.size() > 5)
						mmCifFiles.pop_back();
				}

				auto &pdb_f = *std::get<1>(*ci);

				// Check to see if it is any use to continue with this structure

				auto pdb_chem_comp = pdb_f.data().get("chem_comp");
				if (not pdb_chem_comp)
					continue;

				bool none = true;
				for (const auto &[comp_id] : pdb_chem_comp->rows<std::string>("id"))
				{
					if (ligands[comp_id])
					{
						none = false;
						break;
					}
				}

				if (none)
				{
					// no_compounds.insert(pdb_id);

					if (cif::VERBOSE > 0)
						std::cerr << "This structure does not contain any transplantable compound" << std::endl;

					continue;
				}

				mmcif::Structure pdb_structure(pdb_f);

				// if (not validateHit(pdb_structure, hit))
				// {
				// 	std::cerr << "invalid fasta for hit " << hit.mDefLine << std::endl;
				// 	exit(1);
				// }

				auto af_res = getResiduesForChain(af_structure, "A");
				auto pdb_res = getResiduesForChain(pdb_structure, chain_id);

				if (pdb_res.size() == 0)
				{
					std::cerr << "Missing chain " << chain_id << " in " << pdb_id << std::endl;
					continue;
				}

				for (auto &hsp : hit.mHsps)
				{
					if (cif::VERBOSE > 0)
						std::cerr << "hsp, identity " << std::fixed << std::setprecision(2) << hsp.identity() << " length " << hsp.length() << std::endl;

					if (hsp.length() < minAlignmentLength)
					{
						if (cif::VERBOSE > 0)
							std::cerr << "hsp not long enough" << std::endl;
						continue;
					}

					if (hsp.identity() < minHspIdentity)
					{
						if (cif::VERBOSE > 0)
							std::cerr << "hsp not identical enough" << std::endl;
						continue;
					}

					auto &&[af_ix_trimmed, pdb_ix_trimmed] = getTrimmedIndicesForHsp(hsp);

					// sanity check, happened unfortunately when the fasta was out of sync with the real PDB
					if (pdb_ix_trimmed.back() >= pdb_res.size())
					{
						std::cerr << "Probably incorrect fasta entry for " << hit.mDefLine << std::endl;
						continue;
					}

					assert(af_ix_trimmed.size() == pdb_ix_trimmed.size());

					std::vector<Point> af_ca_trimmed, pdb_ca_trimmed;
					for (size_t i = 0; i < af_ix_trimmed.size(); ++i)
					{
						assert(af_ix_trimmed[i] < af_res.size());
						assert(pdb_ix_trimmed[i] < pdb_res.size());

						auto af_ca = af_res[af_ix_trimmed[i]]->atomByID("CA");
						if (not af_ca)
							continue;

						auto pdb_ca = pdb_res[pdb_ix_trimmed[i]]->atomByID("CA");
						if (not pdb_ca)
							continue;

						af_ca_trimmed.push_back(af_ca.location());
						pdb_ca_trimmed.push_back(pdb_ca.location());
					}

					if (af_ca_trimmed.size() < af_ix_trimmed.size() and cif::VERBOSE > 0)
						std::cerr << "Nr of missing CA: " << (af_ix_trimmed.size() - af_ca_trimmed.size()) << std::endl;
					
					if (af_ca_trimmed.empty())
					{
						if (cif::VERBOSE > 0)
							std::cerr << "No CA atoms mapped, skipping" << std::endl;
						continue;
					}

					double rmsd = Align(af_structure, pdb_structure, af_ca_trimmed, pdb_ca_trimmed);

					json r_hsp{
						{"pdb_id", pdb_id},
						{"pdb_asym_id", pdb_res.front()->asymID()},
						{"identity", hsp.identity()},
						{"alignment_length", hsp.length()},
						{"rmsd", rmsd}};

					for (auto &res : pdb_structure.nonPolymers())
					{
						auto comp_id = res.compoundID();

						Ligand ligand = ligands[comp_id];

						if (not ligand)
							continue;

						// Find the atoms nearby in the AF chain for this residue
						auto &&[pdb_near_r, af_near_r] = selectAtomsNearResidue(
							pdb_res, pdb_ix_trimmed,
							af_res, af_ix_trimmed, res.atoms(), maxLigandBackboneDistance, ligand);

						if (pdb_near_r.size() == 0)
						{
							if (cif::VERBOSE > 0)
								std::cerr << "There are no atoms found near residue " << res << std::endl;
							continue;
						}

						if (cif::VERBOSE > 1)
							std::cerr << "Found " << pdb_near_r.size() << " atoms nearby" << std::endl;

						// realign based on these nearest atoms.
						if (pdb_near_r.size() > 3)
							rmsd = Align(af_structure, pdb_structure, af_near_r, pdb_near_r);
						else if (cif::VERBOSE > 0)
						{
							rmsd = 0;
							if (cif::VERBOSE > 0)
								std::cerr << "There are not enough atoms found near residue " << res << " to fine tune rotation" << std::endl;
						}

						auto analogue = ligand.analogueID();
						if (analogue.empty())
							analogue = comp_id;

						// check to see if the ligand is unique enough
						const auto &[unique, replace_id] = isUniqueLigand(af_structure, minSeparationDistance, res, analogue);

						switch (unique)
						{
							case UniqueType::Seen:
							{
								if (cif::VERBOSE > 0)
									std::cerr << "Residue " << res << " is not unique enough" << std::endl;
								continue;
							}

							case UniqueType::MoreAtoms:
							{
								if (cif::VERBOSE > 0)
								{
									auto &rep_res = af_structure.getResidue(replace_id);
									std::cerr << "Residue " << res << " has more atoms than the first transplant " << rep_res << std::endl;

									af_structure.removeResidue(rep_res);
								}
								break;
							}

							case UniqueType::Unique:
								break;
						}

						// Calculate clash score for new ligand

						std::vector<CAtom> resAtoms;

						for (auto &atom : res.atoms())
						{
							if (ligand.drops(atom.labelAtomID()))
								continue;

							mmcif::AtomTypeTraits att(atom.type());

							int formal_charge = atom.charge();

							if (formal_charge == 0 and att.isMetal() and res.atoms().size() == 1)
							{
								auto compound = mmcif::CompoundFactory::instance().create(comp_id);
								if (compound)
									formal_charge = compound->formalCharge();
							}

							try
							{
								resAtoms.emplace_back(att.type(), atom.location(), formal_charge, atom.labelSeqID(), atom.labelAtomID());
							}
							catch (const std::exception &ex)
							{
								auto compound = mmcif::CompoundFactory::instance().create(att.symbol());
								if (compound)
									formal_charge = compound->formalCharge();
								else
									formal_charge = 0;

								resAtoms.emplace_back(att.type(), atom.location(), formal_charge, atom.labelSeqID(), atom.labelAtomID());
							}
						}

						auto &&[polyAtomCount, clashInfo] = CalculateClashScore(polyAtoms, resAtoms, maxDistance);

						if (polyAtomCount == 0)
						{
							if (cif::VERBOSE > 0)
								std::cerr << "Residue " << res << " skipped because there are no atoms nearby in the polymer" << std::endl;
							continue;
						}

						auto entity_id = af_structure.createNonPolyEntity(comp_id);
						auto asym_id = af_structure.createNonpoly(entity_id, res.atoms());

						r_hsp["transplants"].push_back({
							{"compound_id", comp_id},
							// {"entity_id", entity_id},
							{"asym_id", asym_id},
							{"pdb_asym_id", res.asymID()},
							{"pdb_auth_asym_id", res.authAsymID()},
							{"pdb_auth_seq_id", res.authSeqID()},
							{"rmsd", rmsd},
							{"analogue_id", analogue},
							{"clash", clashInfo}
						});

						if (not res.authInsCode().empty())
							r_hsp["transplants"].back().emplace("pdb_auth_ins_code", res.authInsCode());
						else
							r_hsp["transplants"].back().emplace("pdb_auth_ins_code", nullptr);

						// copy any struct_conn record that might be needed

						auto &pdb_struct_conn = pdb_structure.category("struct_conn");
						auto &af_struct_conn = af_structure.category("struct_conn");

						for (auto atom : res.atoms())
						{
							for (auto conn : pdb_struct_conn.find(
									 ("ptnr1_label_asym_id"_key == atom.labelAsymID() and "ptnr1_label_atom_id"_key == atom.labelAtomID()) or
									 ("ptnr2_label_asym_id"_key == atom.labelAsymID() and "ptnr2_label_atom_id"_key == atom.labelAtomID())))
							{
								std::string a_type, a_comp;
								if (conn["ptnr1_label_asym_id"].as<std::string>() == atom.labelAsymID() and
									conn["ptnr1_label_atom_id"].as<std::string>() == atom.labelAtomID())
								{
									a_type = conn["ptnr2_label_atom_id"].as<std::string>();
									a_comp = conn["ptnr2_label_comp_id"].as<std::string>();
								}
								else
								{
									a_type = conn["ptnr1_label_atom_id"].as<std::string>();
									a_comp = conn["ptnr1_label_comp_id"].as<std::string>();
								}

								// locate the corresponding atom in the af structure
								auto a_a = af_structure.getAtomByPositionAndType(atom.location(), a_type, a_comp);

								if (not a_a)
								{
									if (cif::VERBOSE > 0)
										std::cerr << "Could not create a connection to " << atom << std::endl;
									continue;
								}

								auto conn_type = conn["conn_type_id"].as<std::string>();

								af_struct_conn.emplace({
									{"id", af_struct_conn.getUniqueID(conn_type)},
									{"conn_type_id", conn_type},
									{"ptnr1_label_asym_id", asym_id},
									{"ptnr1_label_comp_id", res.compoundID()},
									{"ptnr1_label_seq_id", "."},
									{"ptnr1_label_atom_id", atom.labelAtomID()},
									{"ptnr1_symmetry", "1_555"},
									{"ptnr2_label_asym_id", a_a.labelAsymID()},
									{"ptnr2_label_comp_id", a_a.labelCompID()},
									{"ptnr2_label_seq_id", a_a.labelSeqID()},
									{"ptnr2_label_atom_id", a_a.labelAtomID()},
									{"ptnr1_auth_asym_id", asym_id},
									{"ptnr1_auth_comp_id", res.compoundID()},
									{"ptnr1_auth_seq_id", "1"},
									{"ptnr1_auth_atom_id", atom.authAtomID()},
									{"ptnr2_auth_asym_id", a_a.labelAsymID()},
									{"ptnr2_auth_comp_id", a_a.labelCompID()},
									{"ptnr2_auth_seq_id", a_a.labelSeqID()},
									{"ptnr2_auth_atom_id", a_a.authAtomID()},
									{"ptnr2_symmetry", "1_555"},
									{"pdbx_dist_value", Distance(a_a, atom)}});
							}
						}

						// now fix up the newly created residue
						ligand.modify(af_structure, asym_id);

						if (cif::VERBOSE > 0)
							std::cerr << "Created asym " << asym_id << " for " << res << std::endl;
					}

					if (not r_hsp["transplants"].empty())
						hits.push_back(std::move(r_hsp));
				}
			}
			catch (const std::exception &e)
			{
				std::throw_with_nested(std::runtime_error("Error when processing " + pdb_id + " for " + afID));
			}
		}
	}

	af_structure.cleanupEmptyCategories();

	af_structure.datablock().add_software("alphafill", "model annotation", kVersionNumber, kBuildDate);

	if (vm.count("output"))
	{
		fs::path output = vm["output"].as<std::string>();

		if (output.has_parent_path() and not fs::exists(output.parent_path()))
			fs::create_directories(output.parent_path());

		f.save(output);

		// if (output.extension() == ".gz")
		// 	output = output.stem();

		output.replace_extension(".json");
		std::ofstream outfile(output);
		outfile << result;
	}
	else
		f.save(std::cout);

	return 0;
}
