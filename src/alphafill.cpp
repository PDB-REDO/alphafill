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

int create_blast_index()
{
	auto &config = mcfp::config::instance();

	fs::path pdbDir = config.get<std::string>("pdb-dir");

	// fs::path ligandsFile = config.get<std::string>("ligands");
	// if (not fs::exists(ligandsFile))
	// 	throw std::runtime_error("Ligands file not found");

	// LigandsTable ligands(ligandsFile);

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
	blocking_queue<std::tuple<std::string,std::string,std::string>> q2;

	std::vector<std::thread> t;
	std::vector<std::vector<std::string>> t_result(nrOfThreads);

	fs::path tmpFile = config.get<std::string>("pdb-fasta") + ".tmp";
	std::ofstream tmpFastA(tmpFile);

	if (not tmpFastA.is_open())
		throw std::runtime_error("Could not open temporary file for writing (" + tmpFile.string() + ")");

	// The writing thread
	std::thread tc([&q2,&tmpFastA]()
	{
		for (;;)
		{
			const auto &[pdb_id, entity_id, seq] = q2.pop();
			if (pdb_id.empty())
				break;
			
			tmpFastA << ">pdb-entity|" << pdb_id << '|' << entity_id << std::endl;

			size_t o = 0;
			for (char l : seq)
			{
				if (std::isspace(l))
					continue;
				
				tmpFastA << l;
				if (++o == 80)
				{
					o = 0;
					tmpFastA << std::endl;
				}
			}

			if (o != 0)
				tmpFastA << std::endl;
		}
	});

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

					auto &entity_poly = db["entity_poly"];

					for (const auto &[entity_id, type, seq] :
						entity_poly.find<std::string,std::string,std::string>(cif::key("type") == "polypeptide(L)",
							"entity_id", "type", "pdbx_seq_one_letter_code_can"))
					{
						q2.push({pdbID, entity_id, seq});
					}
				}
				catch(const std::exception& e)
				{
					std::cerr << e.what() << std::endl;
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
		std::cerr << "Error in walking files: " << ex.what() << std::endl;
		exit(1);
	}

	// signal end
	q1.push({});

	for (auto &ti : t)
		ti.join();

	q2.push({});

	tc.join();

	// replace the old one
	std::error_code ec;
	fs::path fastaFile = config.get<std::string>("pdb-fasta");

	if (fs::exists(fastaFile, ec))
		fs::remove(fastaFile, ec);
	
	if (ec)
		throw std::runtime_error("Could not replace fasta file: " + ec.message());
	
	fs::rename(tmpFile, fastaFile, ec);

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

	std::regex rx(R"(>pdb-entity\|(\w{4,})\|(\w+)( .+)?)");
	if (not std::regex_match(line, rx))
		throw std::runtime_error("The first line in the blast index file does not seem to be correct, please re-create a blast index using the create-blast-index command");
}

zeep::json::element calculatePAEScore(const std::vector<cif::mm::residue *> &af_res, std::vector<CAtom> &atoms, float maxDistance, const PAE_matrix &pae)
{
	auto maxDistanceSq = maxDistance * maxDistance;

	std::vector<size_t> index;

	for (size_t i = 0; i < af_res.size(); ++i)
	{
		for (auto &a : af_res[i]->atoms())
		{
			auto l = a.get_location();
			bool near = false;

			for (auto &b : atoms)
			{
				if (distance_squared(l, b.pt) < maxDistanceSq)
				{
					near = true;
					break;
				}
			}

			if (near)
			{
				index.push_back(i);
				break;
			}
		}
	}

	zeep::json::element result;
	auto &pae_s = result["pae"];

	std::vector<float> s(32);

	for (size_t i = 0; i < index.size(); ++i)
	{
		std::vector<uint8_t> v(index.size());

		for (size_t j = 0; j < index.size(); ++j)
			v[j] = pae(index[i], index[j]);
		pae_s.push_back(v);
	}

	for (size_t i = 0; i < index.size(); ++i)
	{
		for (size_t j = i + 1; j < index.size(); ++j)
		{
			auto v = pae(index[i], index[j]);
			if (v < pae(index[j], index[i]))
				v = pae(index[j], index[i]);

			for (auto k = v; k < 32; ++k)
				s[k] += 1;
		}
	}

	float N = (index.size() * (index.size() - 1)) / 2.0f;
	for (float &v : s)
		v /= N;

	result["scores"] = s;

// std::cout << std::setw(1) << result << std::endl;

	return result;
}

// --------------------------------------------------------------------

zeep::json::element alphafill(cif::datablock &db, const std::vector<PAE_matrix> &v_pae, alphafill_progress_cb &&progress)
{
	using namespace std::literals;
	using namespace cif::literals;

	auto &config = mcfp::config::instance();

	std::string fasta = config.get<std::string>("pdb-fasta");
	fs::path pdbDir = config.get<std::string>("pdb-dir");

	// --------------------------------------------------------------------

	fs::path ligandsFile = config.get<std::string>("ligands");
	if (not fs::exists(ligandsFile))
		throw std::runtime_error("Ligands file not found");

	LigandsTable ligands(ligandsFile);

	float maxDistance = config.get<float>("clash-distance-cutoff");

	// --------------------------------------------------------------------

	cif::iset pdbIDsContainingLigands;
	if (config.has("test-pdb-id"))
		pdbIDsContainingLigands.emplace(config.get<std::string>("test-pdb-id"));
	else if (config.has("pdb-id-list"))
	{
		std::ifstream file(config.get<std::string>("pdb-id-list"));

		if (not file.is_open())
			throw std::runtime_error("Could not open pdb-id-list " + config.get<std::string>("pdb-id-list"));

		std::string line;
		while (std::getline(file, line))
			pdbIDsContainingLigands.insert(line);
	}

	// --------------------------------------------------------------------

	float minHspIdentity = config.get<float>("min-hsp-identity");
	size_t minAlignmentLength = config.get<int>("min-alignment-length");
	float minSeparationDistance = config.get<float>("min-separation-distance");
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

	const std::regex kIDRx(R"(^>pdb-entity\|(\w{4,})\|(\w+)( .*)?)");

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
			std::cerr << "Blasting:" << std::endl
					  << seq << std::endl
					  << std::endl;

		int threads = config.get<int>("threads");
		if (threads == 0)
			threads = std::thread::hardware_concurrency();

		auto result = BlastP(fasta, seq, "BLOSUM62", 3, 10, true, true, 11, 1, config.get<uint32_t>("blast-report-limit"), threads);

		if (cif::VERBOSE > 0)
			std::cerr << "Found " << result.size() << " hits" << std::endl;

		progress.set_max_1(result.size());

		for (auto &hit : result)
		{
			progress.consumed();

			std::smatch m;
			if (not regex_match(hit.mDefLine, m, kIDRx))
			{
				if (cif::VERBOSE > 0)
					std::cerr << "Could not interpret defline from blast hit:" << std::endl
							  << hit.mDefLine << std::endl;
				continue;
			}

			std::string pdb_id = m[1].str();
			std::string entity_id = m[2].str();

			progress.message(pdb_id);

			if (not(pdbIDsContainingLigands.empty() or pdbIDsContainingLigands.count(pdb_id)))
				continue;

			if (cif::VERBOSE > 0)
				std::cerr << "pdb id: " << pdb_id << '\t' << "entity id: " << entity_id << std::endl;

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
						std::cerr << ex.what() << std::endl;
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
						std::cerr << "This structure does not contain any transplantable compound" << std::endl;

					continue;
				}

				cif::mm::structure pdb_structure(pdb_f);

				// if (not validateHit(pdb_structure, hit))
				// {
				// 	std::cerr << "invalid fasta for hit " << hit.mDefLine << std::endl;
				// 	exit(1);
				// }

				for (auto chain_id : get_chain_ids_for_entity_id(pdb_structure.get_datablock(), entity_id))
				{
					auto pdb_res = get_residues_for_chain_id(pdb_structure, chain_id);

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
								std::cerr << "Nr of missing CA: " << (af_ix_trimmed.size() - af_ca_trimmed.size()) << std::endl;

							if (af_ca_trimmed.empty())
							{
								if (cif::VERBOSE > 0)
									std::cerr << "No CA atoms mapped, skipping" << std::endl;
								continue;
							}

							double rmsd = Align(af_structure, pdb_structure, af_ca_trimmed, pdb_ca_trimmed);

							json r_hsp{
								{ "pdb_id", pdb_id },
								{ "pdb_asym_id", pdb_res.front()->get_asym_id() },
								{ "identity", hsp.identity() },
								{ "alignment_length", hsp.length() },
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
											auto &rep_res = af_structure.get_residue(replace_id);
											std::cerr << "Residue " << res << " has more atoms than the first transplant " << rep_res << std::endl;

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
										polyAtoms.emplace_back(atom);
								}

								auto &&[polyAtomCount, clashInfo] = CalculateClashScore(polyAtoms, resAtoms, maxDistance);

								if (polyAtomCount == 0)
								{
									if (cif::VERBOSE > 0)
										std::cerr << "Residue " << res << " skipped because there are no atoms nearby in the polymer" << std::endl;
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
									hsp_t["pae"] = calculatePAEScore(af_res, resAtoms, maxDistance, pae);

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
												std::cerr << "Could not create a connection to " << atom << std::endl;
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

								if (cif::VERBOSE > 0)
									std::cerr << "Created asym " << asym_id << " for " << res << std::endl;
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
	software.emplace({
		{ "pdbx_ordinal", software.size() + 1 }, // TODO: should we check this ordinal number???
		{ "name", "alphafill" },
		{ "version", kVersionNumber },
		{ "date", kBuildDate },
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
	check_blast_index();

	auto &config = mcfp::config::instance();

	if (config.operands().size() < 2)
	{
		std::cout << "Input file not specified" << std::endl;
		exit(1);
	}

	fs::path xyzin = config.operands()[1];

	std::vector<PAE_matrix> v_pae;
	if (config.has("fetch-pae"))
	{
		auto filename = xyzin.filename();

		if (filename.extension() == ".gz")
			filename.replace_extension("");

		if (filename.extension() == ".cif")
			filename.replace_extension("");

		const auto &[type, af_id, chunk, version] = parse_af_id(filename);

		v_pae = data_service::instance().get_pae(af_id, chunk, version);
	}

	cif::file f(xyzin);
	if (f.empty())
		throw std::runtime_error("Empty cif file?");

	if (config.operands().size() >= 3)
	{
		fs::path output = config.operands()[2];

		if (output.has_parent_path() and not fs::exists(output.parent_path()))
			fs::create_directories(output.parent_path());

		json metadata = alphafill(f.front(), v_pae, my_progress{});

		cif::gzio::ofstream xyzout(output);
		f.save(xyzout);

		metadata["file"] = xyzin.string();

		output.replace_extension(".json");
		std::ofstream outfile(output);
		outfile << metadata;
	}
	else
	{
		alphafill(f.front(), v_pae, my_progress{});
		f.save(std::cout);
	}

	return 0;
}
