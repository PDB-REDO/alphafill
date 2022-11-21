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

#include <set>
#include <ostream>

#include <zeep/json/element.hpp>

/// \brief Remove all asymmetric units from the mmCIF file for \a af_id except for the ones in \a requestedAsyms
void stripCifFile(const std::string &af_id, std::set<std::string> requestedAsyms, float identity, std::ostream &os);

/// \brief Using the stripped mmCIF for \a af_id write out an optimized version calculated by Yasara returning statistics
///
/// \param yasara Path to the yasara executable
/// \param af_id The AlphaFold ID to process
/// \param requestedAsyms The list of asym's to include
/// \param identity The identity cut-off
/// \param os The stream to write the output to
/// \result The result is a json object containing the validation statistics before and after running yasara
zeep::json::element optimizeWithYasara(const std::string &af_id, std::set<std::string> requestedAsyms, std::ostream &os);

/// \brief Merge the yasara output in \a yasara_out into an mmCIF structure \a input writing the result to \a os
///
/// \param input Input mmCIF file, alphafill file stripped to one ligand
/// \param yasara_out The Yasara output file containing only atom_site records
/// \param os The stream to write the output to
/// \result The result is a json object containing the validation statistics before and after running yasara
zeep::json::element mergeYasaraOutput(const std::filesystem::path &input, const std::filesystem::path &yasara_out, std::ostream &os);
