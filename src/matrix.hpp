// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// substitution matrix for multiple sequence alignments

#pragma once

#include "cif++/CifUtils.hpp"

// Some predefined matrices

// PAM250 is used by hssp-nt in aligning the sequences

extern const int8_t kMBlosum45[], kMBlosum50[], kMBlosum62[], kMBlosum80[], kMBlosum90[],
	kMPam250[], kMPam30[], kMPam70[];
extern const float kMPam250ScalingFactor, kMPam250MisMatchAverage;

struct MMtrxStats
{
	double lambda, kappa, entropy, alpha, beta;
};

struct MMatrixData
{
	const char *mName;
	int8_t mGapOpen, mGapExtend;
	const int8_t *mMatrix;
	MMtrxStats mGappedStats, mUngappedStats;
};

extern const MMatrixData kMMatrixData[];

// Dayhoff matrix is used for calculating similarity in HSSP

extern const float kDayhoffData[];

// Simple scoring function using the predefined matrices
template <typename T>
inline T score(const T inMatrix[], uint8_t inAA1, uint8_t inAA2)
{
	T result;

	if (inAA1 >= inAA2)
		result = inMatrix[(inAA1 * (inAA1 + 1)) / 2 + inAA2];
	else
		result = inMatrix[(inAA2 * (inAA2 + 1)) / 2 + inAA1];

	return result;
}
