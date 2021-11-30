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

std::vector<compound> data_service::get_compounds() const
{
	pqxx::work tx(db_connection::instance());

	std::vector<compound> compounds;
	for (auto const& [compound_id, analogue_id, s_count_35, t_count_35, s_count_50, t_count_50, s_count_70, t_count_70 ]:
		tx.stream<std::string,std::string,uint32_t, uint32_t,uint32_t, uint32_t,uint32_t, uint32_t>(
		R"(
		 select x.compound_id,
		        x.analogue_id,
				max(s_count_35),
				max(t_count_35),
				max(s_count_50),
				max(t_count_50),
				max(s_count_70),
				max(t_count_70)
		   from (
		 select t.compound_id,
		        t.analogue_id,
				count(distinct h.af_id) as s_count_35,
				count(distinct t.id) as t_count_35,
				0 as s_count_50,
				0 as t_count_50,
				0 as s_count_70,
				0 as t_count_70
		   from af_pdb_hit h
		   join af_transplant t on h.id = t.hit_id
		  where h.identity >= 0.35
		  group by t.compound_id, t.analogue_id
		  union
		 select t.compound_id,
		        t.analogue_id,
				0 as s_count_35,
				0 as t_count_35,
				count(distinct h.af_id) as s_count_50,
				count(distinct t.id) as t_count_50,
				0 as s_count_70,
				0 as t_count_70
		   from af_pdb_hit h
		   join af_transplant t on h.id = t.hit_id
		  where h.identity >= 0.50
		  group by t.compound_id, t.analogue_id
		  union
		 select t.compound_id,
		        t.analogue_id,
				0 as s_count_35,
				0 as t_count_35,
				0 as s_count_50,
				0 as t_count_50,
				count(distinct h.af_id) as s_count_70,
				count(distinct t.id) as t_count_70
		   from af_pdb_hit h
		   join af_transplant t on h.id = t.hit_id
		  where h.identity >= 0.70
		  group by t.compound_id, t.analogue_id
		  ) x
		  group by x.compound_id, x.analogue_id
	      order by x.compound_id asc)"))
	{
		compounds.emplace_back(compound{ compound_id, analogue_id, { s_count_35, s_count_50, s_count_70 }, { t_count_35, t_count_50, t_count_70 } });
	}

	tx.commit();

	return compounds;
}

std::vector<structure> data_service::get_structures(uint32_t page, uint32_t pageSize) const
{
	pqxx::work tx(db_connection::instance());

	const std::regex rx("AF-(.+?)-F1");

	std::vector<structure> structures;
	for (auto const& [structure_id, hit_count_35, transplant_count_35, distinct_35, hit_count_50, transplant_count_50, distinct_50, hit_count_70, transplant_count_70, distinct_70]:
		tx.stream<std::string,uint32_t, uint32_t, uint32_t,uint32_t, uint32_t, uint32_t,uint32_t, uint32_t, uint32_t>(
		R"(select name,
				max(hit_count_35) as max_hit_count_35,
				max(transplant_count_35) as max_transplant_count_35,
				max(distinct_35) as max_distinct_35,
				max(hit_count_50) as max_hit_count_50,
				max(transplant_count_50) as max_transplant_count_50,
				max(distinct_50) as max_distinct_50,
				max(hit_count_70) as max_hit_count_70,
				max(transplant_count_70) as max_transplant_count_70,
				max(distinct_70) as max_distinct_70
			from (
					select s.name as name,
							count(distinct h.id) as hit_count_35,
							count(distinct t.id) as transplant_count_35,
							count(distinct t.analogue_id) as distinct_35,
							0 as hit_count_50,
							0 as transplant_count_50,
							0 as distinct_50,
							0 as hit_count_70,
							0 as transplant_count_70,
							0 as distinct_70
						from af_structure s
						join af_pdb_hit h on s.id = h.af_id
						join af_transplant t on t.hit_id = h.id
						where h.identity >= 0.35
						group by s.name
						union
					select s.name as name,
							0 as hit_count_35,
							0 as transplant_count_35,
							0 as distinct_35,
							count(distinct h.id) as hit_count_50,
							count(distinct t.id) as transplant_count_50,
							count(distinct t.analogue_id) as distinct_50,
							0 as hit_count_70,
							0 as transplant_count_70,
							0 as distinct_70
						from af_structure s
						join af_pdb_hit h on s.id = h.af_id
						join af_transplant t on t.hit_id = h.id
						where h.identity >= 0.5
						group by s.name
						union
					select s.name as name,
							0 as hit_count_35,
							0 as transplant_count_35,
							0 as distinct_35,
							0 as hit_count_50,
							0 as transplant_count_50,
							0 as distinct_50,
							count(distinct h.id) as hit_count_70,
							count(distinct t.id) as transplant_count_70,
							count(distinct t.analogue_id) as distinct_70
						from af_structure s
						join af_pdb_hit h on s.id = h.af_id
						join af_transplant t on t.hit_id = h.id
						where h.identity >= 0.7
						group by s.name
			) x
			group by x.name
			order by max_hit_count_35 desc, name asc
			offset )" + std::to_string(page * pageSize) + R"( rows
			fetch first )" + std::to_string(pageSize) + R"( rows only)"))
	{
		std::smatch m;
		if (std::regex_match(structure_id, m, rx))
			structures.emplace_back(structure{ m[1], { hit_count_35, hit_count_50, hit_count_70 }, { transplant_count_35, transplant_count_50, transplant_count_70 }, { distinct_35, distinct_50, distinct_70 } });
		else
			structures.emplace_back(structure{ structure_id, { hit_count_35, hit_count_50, hit_count_70 }, { transplant_count_35, transplant_count_50, transplant_count_70 }, { distinct_35, distinct_50, distinct_70 } });
	}

	tx.commit();

	return structures;
}

std::vector<structure> data_service::get_structures_for_compound(const std::string &compound, uint32_t page, uint32_t pageSize) const
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
			where t.analogue_id = )" + tx.quote(compound) + R"(
			   or t.compound_id = )" + tx.quote(compound) + R"(
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

uint32_t data_service::count_structures() const
{
	pqxx::work tx(db_connection::instance());

	auto r = tx.exec1(R"(
		  select count(distinct s.id)
			from af_structure s
			right join af_pdb_hit h on s.id = h.af_id
			right join af_transplant t on t.hit_id = h.id
	)");

	tx.commit();

	return r.front().as<uint32_t>();
}

