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

#include "alphafill.hpp"
#include "blast.hpp"
#include "data-service.hpp"
#include "ligands.hpp"
#include "main.hpp"
#include "queue.hpp"
#include "revision.hpp"
#include "utilities.hpp"
#include "validate.hpp"

#include <cif++.hpp>
#include <mcfp/mcfp.hpp>
#include <zeep/json/element.hpp>

#include <chrono>
#include <fstream>
#include <iomanip>

namespace fs = std::filesystem;

using json = zeep::json::element;

// --------------------------------------------------------------------

// The regex to check and/or read FastA files
const std::regex kIDRx(R"(^>pdb-entity\|(\w{4,})\|(\w+)\|([^ |]+)( .*)?)");

// --------------------------------------------------------------------

using cif::point;
using cif::quaternion;

bool validateHit(cif::datablock &db, const BlastHit &hit)
{
	using namespace cif::literals;

	return getSequenceForStrand(db, hit.mDefLine.substr(6)) == hit.mTarget;
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

	return { ixq, ixt };
}

// --------------------------------------------------------------------

enum class UniqueType
{
	Seen,
	Unique,
	MoreAtoms
};

std::tuple<UniqueType, std::string> isUniqueLigand(const cif::mm::structure &structure, float minDistance, const cif::mm::residue &lig, std::string_view id)
{
	std::tuple<UniqueType, std::string> result{ UniqueType::Unique, "" };

	auto minDistanceSq = minDistance * minDistance;

	std::vector<point> pa;

	for (auto &a : lig.atoms())
		pa.push_back(a.get_location());
	auto ca = cif::centroid(pa);

	for (auto &np : structure.non_polymers())
	{
		if (np.get_compound_id() != id)
			continue;
		std::vector<point> pb;

		for (auto &a : np.atoms())
			pb.push_back(a.get_location());
		auto cb = cif::centroid(pb);

		if (distance_squared(ca, cb) < minDistanceSq)
		{
			if (lig.unique_atoms().size() > np.unique_atoms().size())
				result = { UniqueType::MoreAtoms, np.get_asym_id() };
			else
				result = { UniqueType::Seen, np.get_asym_id() };

			break;
		}
	}

	return result;
}

// --------------------------------------------------------------------
//

int create_index(int argc, char *const argv[])
{
	auto &config = load_and_init_config("alphafill create-index [options]",
		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),
		mcfp::make_option<std::string>("pdb-fasta", "The FastA file containing the PDB sequences"),
		mcfp::make_option<int>("threads,t", std::thread::hardware_concurrency(), "Number of threads to use, zero means all available cores"));

	parse_argv(argc, argv, config);

	// --------------------------------------------------------------------

	if (not config.has("pdb-fasta"))
	{
		std::cout << "fasta file not specified\n";
		return 1;
	}

	if (not config.has("pdb-dir"))
	{
		std::cout << "PDB directory not specified\n";
		return 1;
	}

	// --------------------------------------------------------------------

	fs::path pdbDir = config.get("pdb-dir");
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	// --------------------------------------------------------------------

	size_t N = 0;
	for (fs::directory_iterator iter(pdbDir); iter != fs::directory_iterator(); ++iter)
	{
		if (not iter->is_directory())
			continue;

		if (iter->path().filename().string().length() != 2)
			continue;

		++N;
	}

	cif::progress_bar progress(N, "Generating FastA file");

	int nrOfThreads = config.get<int>("threads");
	if (nrOfThreads == 0)
		nrOfThreads = std::thread::hardware_concurrency();

	blocking_queue<fs::path> q1;
	blocking_queue<std::tuple<std::string, std::string>> q2;

	std::vector<std::thread> t;
	std::vector<std::vector<std::string>> t_result(nrOfThreads);

	fs::path tmpFastAFile = config.get("pdb-fasta") + ".tmp";
	std::ofstream tmpFastA(tmpFastAFile);

	if (not tmpFastA.is_open())
		throw std::runtime_error("Could not open temporary file for writing (" + tmpFastAFile.string() + ")");

	// The writing thread
	std::thread tc([&q2, &tmpFastA]()
		{
		for (;;)
		{
			const auto &[id_line, seq] = q2.pop();
			if (id_line.empty())
				break;
			
			tmpFastA << ">pdb-entity|" << id_line << '\n';

			size_t o = 0;
			for (char l : seq)
			{
				if (std::isspace(l))
					continue;
				
				tmpFastA << l;
				if (++o == 80)
				{
					o = 0;
					tmpFastA << '\n';
				}
			}

			if (o != 0)
				tmpFastA << '\n';
		} });

	for (int i = 0; i < nrOfThreads; ++i)
	{
		t.emplace_back([&q1, &q2, ix = i]()
			{
			for (;;)
			{
				fs::path f = q1.pop();

				if (f.empty())	// sentinel
				{
					q1.push({});
					break;
				}

				try
				{
					cif::file pdbFile(f);
					if (pdbFile.empty())
						throw std::runtime_error("invalid cif file " + f.string());

					auto &db = pdbFile.front();
					auto pdbID = db.name();

					if (cif::iequals(pdbID, "XXXX"))
					{
						std::smatch m;
						std::regex rx(R"(([0-9][a-z0-9]{3,7})_final\.cif)", std::regex::icase);

						pdbID = f.string();

						if (std::regex_search(pdbID, m, rx))
						{
							pdbID = m[1];
							std::clog << "\nInvalid PDB-ID in file " << std::quoted(f.string()) << ", using " << std::quoted(pdbID) << " instead\n";
						}
						else
						{
							std::clog << "\nInvalid PDB-ID in file " << std::quoted(f.string()) << ", skipping\n";
							continue;
						}
					}

					std::vector<std::string> ligands;

					auto &chem_comp = db["chem_comp"];
					for (auto &comp_id : chem_comp.rows<std::string>("id"))
					{
						if (cif::compound_factory::instance().is_known_peptide(comp_id) or
							cif::compound_factory::instance().is_known_base(comp_id) or
							comp_id == "HOH")
						{
							continue;
						}

						ligands.emplace_back(comp_id);
					}

					auto &entity_poly = db["entity_poly"];

					for (const auto &[entity_id, type, seq] :
						entity_poly.find<std::string,std::string,std::string>(cif::key("type") == "polypeptide(L)",
							"entity_id", "type", "pdbx_seq_one_letter_code_can"))
					{
						std::ostringstream os;
						os << pdbID << '|' << entity_id << '|';
						bool first = true;
						for (auto &comp_id : ligands)
						{
							if (not std::exchange(first, false))
								os << ';';
							os << comp_id;
						}

						q2.push({os.str(), seq});
					}
				}
				catch(const std::exception& e)
				{
					std::cerr << e.what() << '\n';
				}
			} });
	}

	try
	{
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
				if (not cif::ends_with(name, ".cif.gz"))
					continue;

				q1.push(file);
			}

			// PDB-REDO layout
			for (fs::recursive_directory_iterator fiter(iter->path()); fiter != fs::recursive_directory_iterator(); ++fiter)
			{
				fs::path file = fiter->path();

				std::string name = file.filename().string();
				if (not cif::ends_with(name, "_final.cif"))
					continue;

				q1.push(file);
			}

			progress.consumed(1);
		}
	}
	catch (const std::exception &ex)
	{
		std::cerr << "Error in walking files: " << ex.what() << '\n';
		return 1;
	}

	// signal end
	q1.push({});

	for (auto &ti : t)
		ti.join();

	q2.push({});

	tc.join();

	// replace the old files
	std::error_code ec;
	fs::path fastaFile = config.get("pdb-fasta");

	if (fs::exists(fastaFile, ec))
		fs::remove(fastaFile, ec);

	if (ec)
		throw std::runtime_error("Could not replace fasta file: " + ec.message());

	fs::rename(tmpFastAFile, fastaFile, ec);

	if (ec)
		throw std::runtime_error("Could not rename fasta file: " + ec.message());

	return 0;
}

void check_blast_index()
{
	auto &config = mcfp::config::instance();

	auto indexFile = config.get("pdb-fasta");
	std::ifstream file(indexFile);
	if (not file.is_open())
		throw std::runtime_error("Could not open blast index file (pdb-fasta option)");

	std::string line;
	if (not getline(file, line) or line.empty())
		throw std::runtime_error("blast index file (pdb-fasta option) seems to be empty");

	if (not std::regex_match(line, kIDRx))
		throw std::runtime_error("The first line in the blast index file does not seem to be correct, please re-create a blast index using the create-index command");
}

// --------------------------------------------------------------------

zeep::json::element alphafill(cif::datablock &db, const std::vector<PAE_matrix> &v_pae, alphafill_progress_cb &&progress)
{
	using namespace std::literals;
	using namespace cif::literals;

	auto &config = mcfp::config::instance();

	std::string fasta = config.get("pdb-fasta");
	fs::path pdbDir = config.get("pdb-dir");

	// --------------------------------------------------------------------

	fs::path ligandsFile = config.get("ligands");
	if (not fs::exists(ligandsFile))
		throw std::runtime_error("Ligands file not found");

	LigandsTable ligands(ligandsFile);

	// --------------------------------------------------------------------

	float minHspIdentity = config.get<float>("min-hsp-identity");
	size_t minAlignmentLength = config.get<int>("min-alignment-length");
	float minSeparationDistance = config.get<float>("min-separation-distance");
	float clashDistance = config.get<float>("clash-distance-cutoff");
	float maxLigandBackboneDistance = config.get<float>("max-ligand-to-backbone-distance");

	// --------------------------------------------------------------------

	// This sucks, kinda... The mmcif_af dictionary does not specify
	// all links required to correctly work with libcifpp...
	if (db.get_validator() == nullptr or (db.get_validator()->name() != "mmcif_pdbx.dic" and db.get_validator()->name() != "mmcif_ma.dic"))
		db.set_validator(&cif::validator_factory::instance()["mmcif_pdbx.dic"]);

	cif::mm::structure af_structure(db, 1, cif::mm::StructureOpenOptions::SkipHydrogen);

	if (af_structure.polymers().empty())
		throw std::runtime_error("Structure file does not seem to contain polymers, perhaps pdbx_poly_seq_scheme is missing?");

	// --------------------------------------------------------------------
	// fetch the (single) chain

	std::string afID = db["entry"].front().get<std::string>("id");

	auto v_t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::ostringstream ss;
	ss << std::put_time(std::gmtime(&v_t), "%F");

	json result = {
		{ "id", afID },
		{ "date", ss.str() },
		{ "alphafill_version", kVersionNumber }
	};

	json &hits = result["hits"] = json::array();

	// keep a LRU cache of mmCIF parsed files
	std::list<std::tuple<std::string, std::shared_ptr<cif::file>>> mmCifFiles;

	progress.set_max_0(db["entity_poly"].size());

	auto pae_i = v_pae.begin();
	PAE_matrix empty_pae;

	for (auto r : db["entity_poly"])
	{
		const PAE_matrix &pae = (pae_i == v_pae.end()) ? empty_pae : *pae_i++;

		auto &&[id, seq] = r.get<std::string, std::string>("entity_id", "pdbx_seq_one_letter_code_can");

		// strip all spaces from the sequence, to be able to check length later on
		seq.erase(remove_if(seq.begin(), seq.end(), [](char aa)
					  { return std::isspace(aa); }),
			seq.end());

		if (cif::VERBOSE > 0)
			std::cerr << "Blasting:\n"
					  << seq << '\n'
					  << '\n';

		int threads = config.get<int>("threads");
		if (threads == 0)
			threads = std::thread::hardware_concurrency();

		auto result = BlastP(fasta, seq,
			config.get("blast-matrix"), config.get<int>("blast-word-size"), config.get<double>("blast-expect"),
			not config.has("blast-no-filter"), not config.has("blast-no-gapped"),
			config.get<int>("blast-gap-open"),
			config.get<int>("blast-gap-extend"), config.get<uint32_t>("blast-report-limit"), threads);

		if (cif::VERBOSE > 0)
			std::cerr << "Found " << result.size() << " hits\n";

		progress.set_max_1(result.size());

		for (auto &hit : result)
		{
			progress.consumed();

			std::smatch m;
			if (not regex_match(hit.mDefLine, m, kIDRx))
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Could not interpret defline from blast hit:\n"
							  << hit.mDefLine << '\n';
				continue;
			}

			std::string pdb_id = m[1].str();
			std::string entity_id = m[2].str();
			std::string compound_ids = m[3].str();

			progress.message(pdb_id);

			if (not ligands.contains_any(cif::split(compound_ids, ";")))
				continue;

			if (cif::VERBOSE > 0)
				std::cerr << "pdb id: " << pdb_id << '\t' << "entity id: " << entity_id << '\n';

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
					try
					{
						fs::path pdb_path = pdbFileForID(pdbDir, pdb_id);

						auto cf = std::make_shared<cif::file>(pdb_path.string());

						mmCifFiles.emplace_front(pdb_id, cf);

						// PDB-REDO files don't have the correct audit_conform records, sometimes
						if (cf->get_validator() == nullptr or
							cf->get_validator()->name() == "mmcif_ddl" or
							cf->get_validator()->name() == "mmcif_ddl.dic")
						{
							cf->load_dictionary("mmcif_pdbx");
						}

						ci = mmCifFiles.begin();

						if (mmCifFiles.size() > 5)
							mmCifFiles.pop_back();
					}
					catch (const std::exception &ex)
					{
						std::cerr << ex.what() << '\n';
						continue;
					}
				}

				auto &pdb_f = *std::get<1>(*ci);

				// Check to see if it is any use to continue with this structure

				auto pdb_chem_comp = pdb_f.front().get("chem_comp");
				if (not pdb_chem_comp)
					continue;

				bool none = true;
				for (const auto &comp_id : pdb_chem_comp->rows<std::string>("id"))
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
						std::cerr << "This structure does not contain any transplantable compound\n";

					continue;
				}

				cif::mm::structure pdb_structure(pdb_f);

				// if (not validateHit(pdb_structure, hit))
				// {
				// 	std::cerr << "invalid fasta for hit " << hit.mDefLine << '\n';
				// 	exit(1);
				// }

				for (auto chain_id : get_chain_ids_for_entity_id(pdb_structure.get_datablock(), entity_id))
				{
					auto pdb_res = get_residues_for_chain_id(pdb_structure, chain_id);

					if (pdb_res.size() == 0)
					{
						std::cerr << "Missing chain " << chain_id << " in " << pdb_id << '\n';
						continue;
					}

					for (auto &hsp : hit.mHsps)
					{
						if (cif::VERBOSE > 0)
							std::cerr << "hsp, identity " << std::fixed << std::setprecision(2) << hsp.identity() << " length " << hsp.length() << '\n';

						if (hsp.length() < minAlignmentLength)
						{
							if (cif::VERBOSE > 0)
								std::cerr << "hsp not long enough\n";
							continue;
						}

						if (hsp.identity() < minHspIdentity)
						{
							if (cif::VERBOSE > 0)
								std::cerr << "hsp not identical enough\n";
							continue;
						}

						auto &&[af_ix_trimmed, pdb_ix_trimmed] = getTrimmedIndicesForHsp(hsp);

						// sanity check, happened unfortunately when the fasta was out of sync with the real PDB
						if (pdb_ix_trimmed.back() >= pdb_res.size())
						{
							std::cerr << "Probably incorrect fasta entry for " << hit.mDefLine << '\n';
							continue;
						}

						assert(af_ix_trimmed.size() == pdb_ix_trimmed.size());

						// Loop over each asym containing this entity poly
						for (const auto &af_asym_id : db["struct_asym"].find<std::string>("entity_id"_key == id, "id"))
						{
							auto af_res = get_residuesForAsymID(af_structure, af_asym_id);
							if (af_res.size() != seq.length())
								throw std::runtime_error("Something is wrong with the input file, the number of residues for chain A is not equal to the number in pdbx_seq_one_letter_code_can");

							if (not v_pae.empty() and pae.dim_m() != af_res.size())
								throw std::runtime_error("The supplied PAE data is inconsistent with the residues in the AlphaFold structure for asym ID " + af_asym_id);

							std::vector<point> af_ca_trimmed, pdb_ca_trimmed;
							for (size_t i = 0; i < af_ix_trimmed.size(); ++i)
							{
								assert(af_ix_trimmed[i] < af_res.size());
								assert(pdb_ix_trimmed[i] < pdb_res.size());

								if (af_ix_trimmed[i] >= af_res.size())
									throw std::runtime_error("Mismatch between pdbx_seq_one_letter_code_can and actual residues in input file");
								if (af_ix_trimmed[i] >= af_res.size() or pdb_ix_trimmed[i] >= pdb_res.size())
									throw std::runtime_error("Something is wrong with the data");

								auto af_ca = af_res[af_ix_trimmed[i]]->get_atom_by_atom_id("CA");
								if (not af_ca)
									continue;

								auto pdb_ca = pdb_res[pdb_ix_trimmed[i]]->get_atom_by_atom_id("CA");
								if (not pdb_ca)
									continue;

								af_ca_trimmed.push_back(af_ca.get_location());
								pdb_ca_trimmed.push_back(pdb_ca.get_location());
							}

							if (af_ca_trimmed.size() < af_ix_trimmed.size() and cif::VERBOSE > 0)
								std::cerr << "Nr of missing CA: " << (af_ix_trimmed.size() - af_ca_trimmed.size()) << '\n';

							if (af_ca_trimmed.empty())
							{
								if (cif::VERBOSE > 0)
									std::cerr << "No CA atoms mapped, skipping\n";
								continue;
							}

							double rmsd = Align(af_structure, pdb_structure, af_ca_trimmed, pdb_ca_trimmed);

							json r_hsp{
								{ "pdb_id", pdb_id },
								{ "pdb_asym_id", pdb_res.front()->get_asym_id() },
								{
									"alignment", {
										{ "length", hsp.length() },
										{ "af-begin", hsp.mQueryStart },
										{ "pdb-begin", hsp.mTargetStart },
										{ "identity", hsp.identity() }
									}
								},
								{ "rmsd", rmsd }
							};

							for (auto &res : pdb_structure.non_polymers())
							{
								auto comp_id = res.get_compound_id();

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
										std::cerr << "There are no atoms found near residue " << res << '\n';
									continue;
								}

								if (cif::VERBOSE > 1)
									std::cerr << "Found " << pdb_near_r.size() << " atoms nearby\n";

								// realign based on these nearest atoms.
								if (pdb_near_r.size() > 3)
									rmsd = Align(af_structure, pdb_structure, af_near_r, pdb_near_r);
								else if (cif::VERBOSE > 0)
								{
									rmsd = 0;
									if (cif::VERBOSE > 0)
										std::cerr << "There are not enough atoms found near residue " << res << " to fine tune rotation\n";
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
											std::cerr << "Residue " << res << " is not unique enough\n";
										continue;
									}

									case UniqueType::MoreAtoms:
									{
										if (cif::VERBOSE > 0)
										{
											auto &rep_res = af_structure.get_residue(replace_id);
											std::cerr << "Residue " << res << " has more atoms than the first transplant " << rep_res << '\n';

											af_structure.remove_residue(rep_res);
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
									if (ligand.drops(atom.get_label_atom_id()))
										continue;

									cif::atom_type_traits att(atom.get_type());

									int formal_charge = atom.get_charge();

									if (formal_charge == 0 and att.is_metal() and res.atoms().size() == 1)
									{
										auto compound = cif::compound_factory::instance().create(comp_id);
										if (compound)
											formal_charge = compound->formal_charge();
									}

									try
									{
										resAtoms.emplace_back(att.type(), atom.get_location(), formal_charge, atom.get_label_seq_id(), atom.get_label_atom_id());
									}
									catch (const std::exception &ex)
									{
										auto compound = cif::compound_factory::instance().create(att.symbol());
										if (compound)
											formal_charge = compound->formal_charge();
										else
											formal_charge = 0;

										resAtoms.emplace_back(att.type(), atom.get_location(), formal_charge, atom.get_label_seq_id(), atom.get_label_atom_id());
									}
								}

								// --------------------------------------------------------------------
								// atoms for the clash score calculation
								std::vector<CAtom> polyAtoms;

								for (auto res : af_res)
								{
									for (auto &atom : res->atoms())
									{
										try
										{
											polyAtoms.emplace_back(atom);
										}
										catch (const std::exception &ex)
										{
											polyAtoms.emplace_back(atom.get_type(), atom.get_location(), 0, atom.get_label_seq_id(), atom.get_label_atom_id());
										}
									}
								}

								auto &&[polyAtomCount, clashInfo] = CalculateClashScore(polyAtoms, resAtoms, clashDistance);

								if (polyAtomCount == 0)
								{
									if (cif::VERBOSE > 0)
										std::cerr << "Residue " << res << " skipped because there are no atoms nearby in the polymer\n";
									continue;
								}

								auto entity_id = af_structure.create_non_poly_entity(comp_id);
								auto asym_id = af_structure.create_non_poly(entity_id, res.atoms());

								auto &hsp_t = r_hsp["transplants"].emplace_back(json{
									{ "compound_id", comp_id },
									// {"entity_id", entity_id},
									{ "asym_id", asym_id },
									{ "pdb_asym_id", res.get_asym_id() },
									{ "pdb_auth_asym_id", res.get_auth_asym_id() },
									{ "pdb_auth_seq_id", res.get_auth_seq_id() },
									{ "rmsd", rmsd },
									{ "analogue_id", analogue },
									{ "clash", clashInfo } });

								if (not res.get_pdb_ins_code().empty())
									hsp_t.emplace("pdb_auth_ins_code", res.get_pdb_ins_code());
								else
									hsp_t.emplace("pdb_auth_ins_code", nullptr);

								// Calculate PAE matrix and score for the 'nearby' residues

								if (not pae.empty())
									hsp_t["pae"] = calculatePAEScore(af_res, resAtoms, clashDistance, pae);

								// copy any struct_conn record that might be needed

								auto &pdb_struct_conn = pdb_structure.get_category("struct_conn");
								auto &af_struct_conn = af_structure.get_category("struct_conn");

								for (auto atom : res.atoms())
								{
									for (auto conn : pdb_struct_conn.find(
											 ("ptnr1_label_asym_id"_key == atom.get_label_asym_id() and "ptnr1_label_atom_id"_key == atom.get_label_atom_id()) or
											 ("ptnr2_label_asym_id"_key == atom.get_label_asym_id() and "ptnr2_label_atom_id"_key == atom.get_label_atom_id())))
									{
										std::string a_type, a_comp;
										if (conn["ptnr1_label_asym_id"].as<std::string>() == atom.get_label_asym_id() and
											conn["ptnr1_label_atom_id"].as<std::string>() == atom.get_label_atom_id())
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
										auto a_a = af_structure.get_atom_by_position_and_type(atom.get_location(), a_type, a_comp);

										if (not a_a)
										{
											if (cif::VERBOSE > 0)
												std::cerr << "Could not create a connection to " << atom << '\n';
											continue;
										}

										auto conn_type = conn["conn_type_id"].as<std::string>();

										af_struct_conn.emplace({ { "id", af_struct_conn.get_unique_id(conn_type) },
											{ "conn_type_id", conn_type },
											{ "ptnr1_label_asym_id", asym_id },
											{ "ptnr1_label_comp_id", res.get_compound_id() },
											{ "ptnr1_label_seq_id", "." },
											{ "ptnr1_label_atom_id", atom.get_label_atom_id() },
											{ "ptnr1_label_alt_id", atom.get_label_alt_id() },
											{ "ptnr1_symmetry", "1_555" },
											{ "ptnr2_label_asym_id", a_a.get_label_asym_id() },
											{ "ptnr2_label_comp_id", a_a.get_label_comp_id() },
											{ "ptnr2_label_seq_id", a_a.get_label_seq_id() },
											{ "ptnr2_label_atom_id", a_a.get_label_atom_id() },
											{ "ptnr2_label_alt_id", a_a.get_label_alt_id() },
											{ "ptnr1_auth_asym_id", asym_id },
											{ "ptnr1_auth_comp_id", res.get_compound_id() },
											{ "ptnr1_auth_seq_id", "1" },
											{ "ptnr1_auth_atom_id", atom.get_auth_atom_id() },
											{ "ptnr2_auth_asym_id", a_a.get_label_asym_id() },
											{ "ptnr2_auth_comp_id", a_a.get_label_comp_id() },
											{ "ptnr2_auth_seq_id", a_a.get_label_seq_id() },
											{ "ptnr2_auth_atom_id", a_a.get_auth_atom_id() },
											{ "ptnr2_symmetry", "1_555" },
											{ "pdbx_dist_value", distance(a_a, atom) } });
									}
								}

								// This might not be an alphafold entry, check to see if there's a pdbx_struct_assembly_gen
								if (db.get("pdbx_struct_assembly_gen") != nullptr)
								{
									for (auto r : db["pdbx_struct_assembly_gen"])
									{
										auto asym_id_list = cif::split<std::string>(r["asym_id_list"].as<std::string>(), ",", true);
										if (find(asym_id_list.begin(), asym_id_list.end(), af_res.front()->get_asym_id()) == asym_id_list.end())
											continue;

										asym_id_list.push_back(asym_id);
										r["asym_id_list"] = cif::join(asym_id_list, ",");
										break;
									}
								}

								// now fix up the newly created residue
								ligand.modify(af_structure, asym_id);

								// validation info?
								if (hsp.identity() == 1)
									hsp_t["validation"] = calculateValidationScores(db, pdb_res, af_ix_trimmed, pdb_ix_trimmed,
										af_structure.get_residue(asym_id), res, config.get<float>("max-ligand-to-polymer-atom-distance"), ligand);

								if (cif::VERBOSE > 0)
									std::cerr << "Created asym " << asym_id << " for " << res << '\n';
							}

							if (not r_hsp["transplants"].empty())
								hits.push_back(std::move(r_hsp));
						}
					}
				}
			}
			catch (const std::exception &e)
			{
				std::throw_with_nested(std::runtime_error("Error when processing " + pdb_id + " for " + afID));
			}
		}
	}

	af_structure.cleanup_empty_categories();

	auto &software = af_structure.get_category("software");
	software.emplace({ { "pdbx_ordinal", software.size() + 1 }, // TODO: should we check this ordinal number???
		{ "name", "alphafill" },
		{ "version", kVersionNumber },
		{ "date", kRevisionDate },
		{ "classification", "model annotation" } });

	return result;
}

// --------------------------------------------------------------------

struct my_progress : public alphafill_progress_cb
{
	void set_max_0(size_t in_max) override
	{
	}

	void set_max_1(size_t in_max) override
	{
		if (cif::VERBOSE < 1)
			m_progress.reset(new cif::progress_bar(in_max + 1, "matching"));
	}

	void consumed(size_t n) override
	{
		if (m_progress)
			m_progress->consumed(n);
	}

	void message(const std::string &msg) override
	{
		if (m_progress)
			m_progress->message(msg);
	}

	std::unique_ptr<cif::progress_bar> m_progress;
};

int alphafill_main(int argc, char *const argv[])
{
	auto &config = load_and_init_config("usage: alphafill process [options] <inputfile> [<outputfile>]",

		mcfp::make_option<std::string>("pae-file", "Specify a specific file containing PAE information, default is to use a filename based on inputfile"),

		mcfp::make_option<std::string>("pdb-dir", "Directory containing the mmCIF files for the PDB"),
		mcfp::make_option<std::string>("pdb-fasta", "The FastA file containing the PDB sequences"),

		mcfp::make_option<std::string>("ligands", "af-ligands.cif", "File in CIF format describing the ligands and their modifications"),

		mcfp::make_option<float>("max-ligand-to-backbone-distance", 6, "The max distance to use to find neighbouring backbone atoms for the ligand in the AF structure"),
		mcfp::make_option<float>("min-hsp-identity", 0.25, "The minimal identity for a high scoring pair (note, value between 0 and 1)"),
		mcfp::make_option<int>("min-alignment-length", 85, "The minimal length of an alignment"),
		mcfp::make_option<float>("min-separation-distance", 3.5, "The centroids of two identical ligands should be at least this far apart to count as separate occurrences"),

		mcfp::make_option<float>("clash-distance-cutoff", 4, "The max distance between polymer atoms and ligand atoms used in calculating clash scores"),

		mcfp::make_option<float>("max-ligand-to-polymer-atom-distance", 6,
			"The max distance to use to find neighbouring polymer atoms for the ligand in the AF structure (for validation)"),

		mcfp::make_option<uint32_t>("blast-report-limit", 250, "Number of blast hits to use"),

		mcfp::make_hidden_option<std::string>("blast-matrix", "BLOSUM62", "Blast matrix to use"),
		mcfp::make_hidden_option<int>("blast-word-size", 3, "Blast word size"),
		mcfp::make_hidden_option<double>("blast-expect", 10, "Blast expect cut off"),
		mcfp::make_hidden_option("blast-no-filter", "Blast option for filter, default is to use low complexity filter"),
		mcfp::make_hidden_option("blast-no-gapped", "Blast option for gapped, default is to do gapped"),
		mcfp::make_hidden_option<int>("blast-gap-open", 11, "Blast penalty for gap open"),
		mcfp::make_hidden_option<int>("blast-gap-extend", 1, "Blast penalty for gap extend"),

		mcfp::make_option<int>("threads,t", std::thread::hardware_concurrency(), "Number of threads to use, zero means all available cores"),

		mcfp::make_hidden_option<std::string>("custom-dir", (fs::temp_directory_path() / "alphafill").string(), "Directory for custom built entries")

	);

	parse_argv(argc, argv, config);

	if (config.operands().size() > 2 or config.operands().size() < 1)
	{
		std::cerr << config << '\n';
		return 1;
	}

	// --------------------------------------------------------------------

	if (not config.has("pdb-fasta"))
	{
		std::cout << "fasta file not specified\n";
		return 1;
	}

	check_blast_index();

	if (not config.has("pdb-dir"))
	{
		std::cout << "PDB directory not specified\n";
		return 1;
	}

	// --------------------------------------------------------------------

	fs::path xyzin = config.operands().front();

	cif::file f(xyzin);
	if (f.empty())
	{
		std::cerr << "Empty cif file?\n";
		return 1;
	}

	fs::path paein;

	if (config.has("pae-file"))
		paein = config.get("pae-file");
	else
	{
		auto filename = xyzin.filename();

		if (filename.extension() == ".gz")
			filename.replace_extension();

		if (filename.extension() == ".cif")
			filename.replace_extension();

		const auto &[type, af_id, chunk, version] = parse_af_id(filename);

		// paein = xyzin.parent_path() / std::format("AF-{}-F{}-predicted_aligned_error_v{}.json", af_id, chunk, version);
		paein = xyzin.parent_path() / cif::format("AF-%s-F%d-predicted_aligned_error_v%d.json", af_id, chunk, version).str();
	}

	std::vector<PAE_matrix> v_pae;
	if (fs::exists(paein))
		v_pae = data_service::load_pae_from_file(paein);

	json metadata = alphafill(f.front(), v_pae, my_progress{});

	if (config.operands().size() == 2)
	{
		fs::path output = config.operands().back();

		if (output.has_parent_path() and not fs::exists(output.parent_path()))
			fs::create_directories(output.parent_path());

		cif::gzio::ofstream xyzout(output);
		f.save(xyzout);

		metadata["file"] = xyzin.string();

		output.replace_extension(".json");
		std::ofstream outfile(output);
		outfile << metadata;
	}
	else
		f.save(std::cout);

	return 0;
}
