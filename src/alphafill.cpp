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

#include <cif++/Structure.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "blast.hpp"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;

// --------------------------------------------------------------------

namespace
{
std::string gVersionNr, gVersionDate;
}

void load_version_info()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)"),
		rxVersionNr2(R"(mkdssp-version: (\d+(?:\.\d+)+))");

#include "revision.hpp"

	struct membuf : public std::streambuf
	{
		membuf(char *data, size_t length) { this->setg(data, data, data + length); }
	} buffer(const_cast<char *>(kRevision), sizeof(kRevision));

	std::istream is(&buffer);

	std::string line;

	while (getline(is, line))
	{
		std::smatch m;

		if (std::regex_match(line, m, rxVersionNr))
		{
			gVersionNr = m[1];
			if (m[2].matched)
				gVersionNr += '*';
			continue;
		}

		if (std::regex_match(line, m, rxVersionDate))
		{
			gVersionDate = m[1];
			continue;
		}

		// always the first, replace with more specific if followed by the other info
		if (std::regex_match(line, m, rxVersionNr2))
		{
			gVersionNr = m[1];
			continue;
		}
	}
}

std::string get_version_nr()
{
	return gVersionNr /* + '/' + cif::get_version_nr()*/;
}

std::string get_version_date()
{
	return gVersionDate;
}

std::string get_version_string()
{
	return gVersionNr + " " + gVersionDate;
}

// --------------------------------------------------------------------

using mmcif::Point;
using mmcif::Quaternion;

// std::vector<Point> getCAlphaForChain(cif::Datablock &db, const std::string &asym_id)
// {
// 	using namespace cif::literals;

// 	std::vector<Point> result;

// 	auto &atoms = db["atom_site"];
// 	for (const auto &[x, y, z] : atoms.find<float, float, float>("auth_asym_id"_key == asym_id and "label_atom_id"_key == "CA", "Cartn_x", "Cartn_y", "Cartn_z"))
// 		result.emplace_back(x, y, z);

// 	return result;
// }

std::vector<mmcif::Residue*> getResiduesForChain(mmcif::Structure &structure, const std::string &asym_id)
{
	std::vector<mmcif::Residue*> result;

	for (auto &poly: structure.polymers())
	{
		if (poly.asymID() != asym_id)
			continue;
		
		for (auto &res: poly)
			result.emplace_back(&res);
	}

	return result;
}

std::vector<Point> getCAlphaForChain(const std::vector<mmcif::Residue*> &residues)
{
	std::vector<Point> result;

	for (auto res: residues)
		result.push_back(res->atomByID("CA").location());

	return result;
}

std::tuple<std::vector<size_t>, std::vector<size_t>> getTrimmedIndicesForHsp(const BlastHsp &hsp)
{
	std::vector<size_t> ixq, ixt;

	assert(hsp.mAlignedQuery.length() == hsp.mAlignedTarget.length());

	size_t qix = hsp.mQueryStart;
	size_t tix = hsp.mTargetStart;

	auto &qa = hsp.mAlignedQuery;
	auto &ta = hsp.mAlignedTarget;

	for (size_t i = 0; i < hsp.mAlignedQuery.length(); ++i)
	{
		if (is_gap(qa[qix]))
		{
			++tix;
			continue;
		}

		if (is_gap(ta[tix]))
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

double CalculateRMSD(const std::vector<mmcif::Point> &pa, const std::vector<mmcif::Point> &pb)
{
	return RMSd(pa, pb);
}

void Align(mmcif::Structure &a, mmcif::Structure &b,
	std::vector<Point> &cAlphaA, std::vector<Point> &cAlphaB)
{
	auto ta = Centroid(cAlphaA);

	if (cif::VERBOSE)
		std::cerr << "translate A: " << -ta << std::endl;

	for (auto &pt : cAlphaA)
		pt -= ta;

	auto tb = Centroid(cAlphaB);

	if (cif::VERBOSE)
		std::cerr << "translate B: " << -tb << std::endl;

	for (auto &pt : cAlphaB)
		pt -= tb;

	auto rotation = AlignPoints(cAlphaA, cAlphaB);
	// for (auto &pt : cAlphaB)
	// 	pt.rotate(rotation);

	if (cif::VERBOSE)
	{
		const auto & [ angle, axis ] = mmcif::QuaternionToAngleAxis(rotation);
		std::cerr << "rotation: " << angle << " degrees rotation around axis " << axis << std::endl;
	}

	a.translate(-ta);
	b.translate(-tb);
	b.rotate(rotation);

	if (cif::VERBOSE)
		std::cerr << "RMSd: " << CalculateRMSD(cAlphaA, cAlphaB) << std::endl;
}

std::tuple<std::vector<Point>, std::vector<Point>> selectAtomsNearResidue(
	const std::vector<mmcif::Residue*> &pdb, const std::vector<size_t> &pdb_ix,
	const std::vector<mmcif::Residue*> &af, const std::vector<size_t> &af_ix,
	const std::vector<mmcif::Atom> &residue, float maxDistance)
{
	std::vector<Point> ra, rb;

	assert(pdb_ix.size() == af_ix.size());

	for (size_t i = 0; i < pdb_ix.size(); ++i)
	{
		bool nearby = false;

		for (const char *atom_id: { "C", "CA", "N", "O" })
		{
			try
			{
				auto atom = pdb[pdb_ix[i]]->atomByID(atom_id);
				for (auto &b: residue)
				{
					if (Distance(atom, b) <= maxDistance)
					{
						nearby = true;
						break;
					}
				}
			}
			catch(const std::exception& e)
			{
			}
			
			if (nearby)
				break;
		}

		if (not nearby)
			continue;

		for (const char *atom_id: { "C", "CA", "N", "O" })
		{
			try
			{
				auto pt_a = pdb[pdb_ix[i]]->atomByID(atom_id).location();
				auto pt_b = af[af_ix[i]]->atomByID(atom_id).location();

				ra.push_back(pt_a);
				rb.push_back(pt_b);
			}
			catch(const std::exception& e)
			{
			}
		}
	}

	return { ra, rb };
}

// --------------------------------------------------------------------

int a_main(int argc, const char *argv[])
{
	using namespace std::literals;
	using namespace cif::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);
	visible_options.add_options()
		("pdb-fasta",		po::value<std::string>(),	"The FastA file containing the PDB sequences")
		("pdb-dir",			po::value<std::string>(),	"Directory containing the mmCIF files for the PDB")
		("transplantable",	po::value<std::string>(),	"Semicolon separated list of transplantable residues")
		("compounds",		po::value<std::string>(),	"Location of the components.cif file from CCD")
		("components",		po::value<std::string>(),	"Location of the components.cif file from CCD, alias")
		("extra-compounds",	po::value<std::string>(),	"File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file")
		("mmcif-dictionary",po::value<std::string>(),	"Path to the mmcif_pdbx.dic file to use instead of default")
		("config",			po::value<std::string>(),	"Config file")
		("help,h",										"Display help message")
		("version",										"Print version")
		("verbose,v",									"Verbose output");

	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("xyzin,i",			po::value<std::string>(),	"coordinates file")
		("output,o",		po::value<std::string>(),	"Output to this file")
		("debug,d",			po::value<int>(),			"Debug level (for even more verbose output)")

		("test",										"Run test code")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	fs::path configFile = "alphafill.conf";
	if (vm.count("config"))
		configFile = vm["config"].as<std::string>();

	if (not fs::exists(configFile) and getenv("HOME") != nullptr)
		configFile = fs::path(getenv("HOME")) / ".config" / configFile;
	
	if (fs::exists(configFile))
	{
		std::ifstream cfgFile(configFile);
		if (cfgFile.is_open())
			po::store(po::parse_config_file(cfgFile, visible_options), vm);
	}

	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << get_version_string() << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cout << visible_options << std::endl;
		exit(0);
	}

	if (vm.count("xyzin") == 0)
	{
		std::cout << "Input file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("pdb-fasta") == 0)
	{
		std::cout << "fasta file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("transplantable") == 0)
	{
		std::cout << "List of transplantable residues not specified" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// --------------------------------------------------------------------

	// Load extra CCD definitions, if any

	if (vm.count("compounds"))
		cif::addFileResource("components.cif", vm["compounds"].as<std::string>());
	else if (vm.count("components"))
		cif::addFileResource("components.cif", vm["components"].as<std::string>());

	if (vm.count("extra-compounds"))
		mmcif::CompoundFactory::instance().pushDictionary(vm["extra-compounds"].as<std::string>());

	// And perhaps a private mmcif_pdbx dictionary

	if (vm.count("mmcif-dictionary"))
		cif::addFileResource("mmcif_pdbx_v50.dic", vm["mmcif-dictionary"].as<std::string>());

	// --------------------------------------------------------------------

	std::string fasta = vm["pdb-fasta"].as<std::string>();
	if (not fs::exists(fasta))
		throw std::runtime_error("PDB-Fasta file does not exist (" + fasta + ")");

	fs::path pdbDir = vm["pdb-dir"].as<std::string>();
	if (not fs::is_directory(pdbDir))
		throw std::runtime_error("PDB directory does not exist");

	std::set<std::string> transplantable;
	ba::split(transplantable, vm["transplantable"].as<std::string>(), ba::is_any_of(";, "));

	// --------------------------------------------------------------------

	mmcif::File f(vm["xyzin"].as<std::string>());
	mmcif::Structure af_structure(f, 1, mmcif::StructureOpenOptions::SkipHydrogen);

	// if (vm.count("test"))
	// {
	// 	Quaternion q{ 0.0235969, 0.24543, 0.800531, 0.546221 };

	// 	af_structure.rotate(q);
	// 	af_structure.translate({ 1, 2, 3 });

	// 	if (vm.count("output"))
	// 		f.file().save(vm["output"].as<std::string>());
	// 	else
	// 		f.file().save(std::cout);
	// 	return 0;
	// }

	// fetch the (single) chain

	const std::regex kIDRx(R"(^>(\w{4,8})_(\w)( .*)?)");

	auto &db = f.data();

	bool done = false;

	for (auto r : db["entity_poly"])
	{
		auto &&[id, seq] = r.get<std::string, std::string>({"entity_id", "pdbx_seq_one_letter_code"});

		std::cout << "Blasting:" << std::endl
				  << seq << std::endl
				  << std::endl;

		auto result = BlastP(fasta, seq);

		std::cout << "Found " << result.size() << " hits" << std::endl;

		for (auto &hit : result)
		{
			std::cout << hit.mDefLine << '\t'
					  << decode(hit.mTarget) << std::endl;

			std::smatch m;
			if (not regex_match(hit.mDefLine, m, kIDRx))
				continue;

			std::string pdb_id = m[1].str();
			std::string chain_id = m[2].str();

			std::cout << "pdb id: " << pdb_id << '\t' << "chain id: " << chain_id << std::endl;

			try
			{
				mmcif::File pdb_f((pdbDir / pdb_id.substr(1, 2) / (pdb_id + ".cif.gz")).string());
				mmcif::Structure pdb_structure(pdb_f);

				auto af_res = getResiduesForChain(af_structure, "A");
				auto pdb_res = getResiduesForChain(pdb_structure, chain_id);

				if (pdb_res.size() == 0)
				{
					std::cerr << "Missing chain " << chain_id << std::endl;
					continue;
				}

				auto &&[af_ix_trimmed, pdb_ix_trimmed] = getTrimmedIndicesForHsp(hit.mHsps.front());

				std::vector<Point> af_ca_trimmed, pdb_ca_trimmed;
				for (size_t i = 0; i < af_ix_trimmed.size(); ++i)
				{
					try
					{
						auto af_ca = af_res[af_ix_trimmed[i]]->atomByID("CA").location();
						auto pdb_ca = pdb_res[pdb_ix_trimmed[i]]->atomByID("CA").location();

						af_ca_trimmed.push_back(af_ca);
						pdb_ca_trimmed.push_back(pdb_ca);
					}
					catch(const std::exception& e)
					{
					}
				}

				Align(af_structure, pdb_structure, af_ca_trimmed, pdb_ca_trimmed);

				for (auto &np : pdb_structure.nonPolymers())
				{
					auto comp_id = np.compoundID();

					if (not transplantable.count(comp_id))
						continue;

					// fetch the non poly atoms as cif::Rows
					auto &res = pdb_structure.getResidue(np.asymID(), comp_id);

					// Find the atoms nearby in the AF chain for this residue
					auto && [ pdb_near_r, af_near_r ] = selectAtomsNearResidue(
						pdb_res, pdb_ix_trimmed,
						af_res, af_ix_trimmed, res.atoms(), 6.f);

					if (pdb_near_r.size() == 0)
					{
						if (cif::VERBOSE)
							std::cerr << "There are no atoms found near residue " << res << std::endl;
						continue;
					}

					// realign based on these nearest atoms.
					if (pdb_near_r.size() > 3)
						Align(af_structure, pdb_structure, af_near_r, pdb_near_r);
					else if (cif::VERBOSE)
						std::cerr << "There are not enough atoms found near residue " << res << " to fine tune rotation" << std::endl;

					auto entity_id = af_structure.createNonPolyEntity(comp_id);
					auto asym_id = af_structure.createNonpoly(entity_id, res.atoms());

					if (cif::VERBOSE)
						std::cerr << "Created asym " << asym_id << std::endl;

					done = true;
				}
			}
			catch (const std::exception &e)
			{
				std::cout << e.what() << '\n';
			}

			if (done)
				break;
		}
		if (done)
			break;
	}

	if (vm.count("output"))
		f.file().save(vm["output"].as<std::string>());
	else
		f.file().save(std::cout);

	return 0;
}

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cout << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cout << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

int main(int argc, const char *argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::addDataDirectory(DATA_DIR);
#endif
		load_version_info();

		result = a_main(argc, argv);
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
