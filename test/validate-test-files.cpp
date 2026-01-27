#include <mcfp/mcfp.hpp>

#include <cif++.hpp>
#include <zeep/el/object.hpp>

#include <fstream>

int main(int argc, char * const argv[])
{
	using namespace cif::literals;

	int result = 0;

	auto &config = mcfp::config::instance();

	config.init("usage: validate <cif-file> <json-file> --pdb-id=<PDBID> --pdb-asym-id=<PDBASYMID>",
		mcfp::make_option<std::string>("pdb-id", "PDB ID"),
		mcfp::make_option<std::string>("pdb-asym-id", "PDB asym ID")
	);

	std::error_code ec;
	config.parse(argc, argv, ec);

	if (ec or config.operands().size() != 2)
	{
		std::cerr << config << '\n';
		exit(1);
	}

	auto pdb_id = config.get("pdb-id");
	auto pdb_asym_id = config.get("pdb-asym-id");

	// JSON first
	std::ifstream json_in(config.operands().back());
	if (not json_in.is_open())
	{
		std::cerr << "Could not open json file " << std::quoted(config.operands().back()) << '\n';
		exit(1);
	}

	auto jdata = zeep::el::object::parse_JSON(json_in);

	auto hits = jdata["hits"];
	if (not hits.is_array() or hits.empty())
	{
		std::cerr << "Unexpected hits array in json\n";
		exit(1);
	}

	auto first_hit = hits[0];
	if (first_hit["pdb_id"] != pdb_id)
	{
		std::cerr << "First hit pdb_id should be " << pdb_id << '\n';
		result = 1;
	}

	if (first_hit["pdb_asym_id"] != pdb_asym_id)
	{
		std::cerr << "First hit pdb_asym_id should be " << pdb_asym_id << '\n';
		result = 1;
	}

	auto transplants = first_hit["transplants"];
	if (not transplants.is_array() or transplants.empty())
	{
		std::cerr << "Unexpected transplants array in json\n";
		exit(1);
	}

	auto first_transplant = transplants[0];

	auto compound_id = first_transplant["compound_id"].get<std::string>();
	auto asym_id = first_transplant["asym_id"].get<std::string>();
	
	cif::file f(config.operands()[0]);
	auto &non_polies = f.front()["pdbx_nonpoly_scheme"];
	if (non_polies.empty())
	{
		std::cerr << "no non-polymers\n";
		result = 1;
	}

	if (not non_polies.contains("asym_id"_key == asym_id and "mon_id"_key == compound_id))
	{
		std::cerr << "The nonpoly scheme does not have a corresponding entry with compound_id/asym_id\n";
		result = 1;
	}

	return result;
}

