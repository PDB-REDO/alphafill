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

#include <regex>

#include "data-service.hpp"
#include "db-connection.hpp"

// --------------------------------------------------------------------

data_service &data_service::instance()
{
	static data_service s_instance;
	return s_instance;
}

std::vector<compound> data_service::get_compounds(float min_identity) const
{
	pqxx::work tx(db_connection::instance());

	std::vector<compound> compounds;
	for (auto const& [compound_id, analogue_id, s_count, t_count ]:
		tx.stream<std::string,std::string,uint32_t, uint32_t>(
		R"(
		 select t.compound_id,
		        t.analogue_id,
				count(distinct h.af_id) as s_count,
				count(distinct t.id) as t_count
		   from af_pdb_hit h
		   join af_transplant t on h.id = t.hit_id
		  where h.identity >= )" + std::to_string(min_identity) + R"(
		  group by t.compound_id, t.analogue_id
	      order by t.compound_id asc)"))
	{
		compounds.emplace_back(compound{ compound_id, analogue_id, s_count, t_count });
	}

	tx.commit();

	return compounds;
}

std::vector<structure> data_service::get_structures(float min_identity, uint32_t page, uint32_t pageSize) const
{
	pqxx::work tx(db_connection::instance());

	const std::regex rx("AF-(.+?)-F1");

	std::vector<structure> structures;
	for (auto const& [structure_id, hit_count, transplant_count, distinct]:
		tx.stream<std::string,uint32_t, uint32_t, uint32_t>(
		R"(select s.name as name,
				count(distinct h.id) as hit_count,
				count(distinct t.id) as transplant_count,
				count(distinct t.analogue_id) as distinct
			from af_structure s
			join af_pdb_hit h on s.id = h.af_id
			join af_transplant t on t.hit_id = h.id
			where h.identity >= )" + std::to_string(min_identity) + R"(
			group by s.name
			order by hit_count desc, s.name asc
			offset )" + std::to_string(page * pageSize) + R"( rows
			fetch first )" + std::to_string(pageSize) + R"( rows only)"))
	{
		std::smatch m;
		if (std::regex_match(structure_id, m, rx))
			structures.emplace_back(structure{ m[1], hit_count, transplant_count, distinct });
		else
			structures.emplace_back(structure{ structure_id, hit_count, transplant_count, distinct });
	}

	tx.commit();

	return structures;
}

std::vector<structure> data_service::get_structures_for_compound(float min_identity, const std::string &compound, uint32_t page, uint32_t pageSize) const
{
	pqxx::work tx(db_connection::instance());

	const std::regex rx("AF-(.+?)-F1");

	std::vector<structure> structures;
	for (auto const& [structure_id, hit_count, transplant_count, distinct]:
		tx.stream<std::string,uint32_t, uint32_t, uint32_t>(
		R"(select s.name,
		          count(distinct h.id) as hit_count,
				  count(distinct t.id) as transplant_count,
				  count(distinct t.analogue_id) as dist_transplant_count
			 from af_structure s
			 join af_pdb_hit h on s.id = h.af_id
			 join af_transplant t on t.hit_id = h.id
			where (t.analogue_id = )" + tx.quote(compound) + R"(
			   or t.compound_id = )" + tx.quote(compound) + R"()
			  and h.identity >= )" + std::to_string(min_identity) + R"(
			group by s.name
			order by hit_count desc
			offset )" + std::to_string(page * pageSize) + R"( rows
			fetch first )" + std::to_string(pageSize) + R"( rows only)"))
	{
		std::smatch m;
		if (std::regex_match(structure_id, m, rx))
			structures.emplace_back(structure{ m[1], hit_count, transplant_count, distinct });
		else
			structures.emplace_back(structure{ structure_id, hit_count, transplant_count, distinct });
	}

	tx.commit();

	return structures;
}

uint32_t data_service::count_structures(float min_identity) const
{
	pqxx::work tx(db_connection::instance());

	auto r = tx.exec1(R"(
		  select count(distinct s.id)
			from af_structure s
			right join af_pdb_hit h on s.id = h.af_id
			right join af_transplant t on t.hit_id = h.id
			where h.identity >= )" + std::to_string(min_identity)
	);

	tx.commit();

	return r.front().as<uint32_t>();
}

uint32_t data_service::count_structures(float min_identity, const std::string &compound) const
{
	pqxx::work tx(db_connection::instance());

	auto r = tx.exec1(R"(
		  select count(distinct s.id)
			from af_structure s
			right join af_pdb_hit h on s.id = h.af_id
			right join af_transplant t on t.hit_id = h.id
			where h.identity >= )" + std::to_string(min_identity) + R"(
			  and (t.compound_id = )" + tx.quote(compound) + " or t.analogue_id = " + tx.quote(compound) + ")"
	);

	tx.commit();

	return r.front().as<uint32_t>();
}

