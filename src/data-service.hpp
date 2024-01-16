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

#pragma once

#include <atomic>
#include <filesystem>
#include <thread>

#include <zeep/nvp.hpp>

#include <cif++.hpp>

#include "queue.hpp"

// --------------------------------------------------------------------

struct compound
{
	std::string id;
	std::string analogue;
	// std::string name;
	uint32_t count_structures;
	uint32_t count_transplants;

	template<typename Archive>
	void serialize(Archive &ar, unsigned long)
	{
		ar & zeep::make_nvp("id", id)
		   & zeep::make_nvp("analogue", analogue)
		   & zeep::make_nvp("structure-count", count_structures)
		   & zeep::make_nvp("transplant-count", count_transplants);
	}
};

struct structure
{
	std::string name;
	uint32_t count_hits;
	uint32_t count_transplants;
	uint32_t distinct_analogues;

	template<typename Archive>
	void serialize(Archive &ar, unsigned long)
	{
		ar & zeep::make_nvp("name", name)
		   & zeep::make_nvp("hit-count", count_hits)
		   & zeep::make_nvp("transplant-count", count_transplants)
		   & zeep::make_nvp("distinct-analogues", distinct_analogues);
	}
};

enum class CustomStatus
{
	Unknown, Queued, Running, Finished, Error
};

struct status_reply
{
	CustomStatus status;
	std::optional<float> progress;
	std::optional<std::string> message;

	template<typename Archive>
	void serialize(Archive &ar, unsigned long)
	{
		ar & zeep::make_nvp("status", status)
		   & zeep::make_nvp("progress", progress)
		   & zeep::make_nvp("message", message);
	}
};

class data_service
{
  public:
	static data_service &instance();

	~data_service();

	void start_queue(size_t nr_of_threads);

	static int rebuild(const std::string &db_user, const std::filesystem::path &db_dir);

	std::vector<compound> get_compounds(float min_identity) const;
	std::vector<structure> get_structures(float min_identity, uint32_t page, uint32_t pageSize) const;
	std::vector<structure> get_structures_for_compound(float min_identity, const std::string &compound, uint32_t page, uint32_t pageSize) const;

	uint32_t count_structures(float min_identity) const;
	uint32_t count_structures(float min_identity, const std::string &compound) const;

	// On demand services

	bool exists_in_afdb(const std::string &id) const;
	std::tuple<std::filesystem::path,std::string,std::string> fetch_from_afdb(const std::string &id) const;

	status_reply get_status(const std::string &id) const;

	void queue(const std::string &data, const std::optional<std::string> pae, const std::string &id);
	std::string queue_af_id(const std::string &id);

  private:

	data_service();

	void run();

	void process_queued(const std::filesystem::path &xyzin, const std::filesystem::path &paein,
		const std::filesystem::path &xyzout, const std::filesystem::path &jsonout);

	std::filesystem::path m_in_dir;
	std::filesystem::path m_out_dir;
	std::filesystem::path m_work_dir;

	// threads processing AF entries
	std::vector<std::thread> m_threads;
	blocking_queue<std::string,100> m_queue;

	std::mutex m_mutex;
	std::string m_running;
	std::atomic<float> m_progress;
};
