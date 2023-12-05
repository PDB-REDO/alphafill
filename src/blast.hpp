// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
//	Simplified Blast algorithm implementation. Works on
//	FastA formatted files containing proteins.

#pragma once

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

// --------------------------------------------------------------------

using sequence = std::basic_string<uint8_t>;

// 22 real letters and 1 dummy (X is the dummy, B and Z are pseudo letters)
extern const char kResidues[]; // = "ACDEFGHIKLMNPQRSTVWYBZX";
extern const uint8_t kResidueNrTable[];

inline constexpr uint8_t ResidueNr(char inAA)
{
	uint8_t result = 23;

	inAA |= 040;
	if (inAA >= 'a' and inAA <= 'z')
		result = kResidueNrTable[inAA - 'a'];

	return result;
}

inline constexpr bool is_gap(char aa)
{
	return aa == ' ' or aa == '.' or aa == '-';
}

inline constexpr bool is_gap(uint8_t aa)
{
	return aa == '-';
}

sequence encode(const std::string &s);
std::string decode(const sequence &s);

// --------------------------------------------------------------------

class blast_exception : public std::runtime_error
{
  public:
	blast_exception(std::string msg)
		: std::runtime_error(msg.c_str())
	{
	}
};

// --------------------------------------------------------------------

struct BlastHsp
{
	uint32_t mScore;
	uint32_t mQueryStart, mQueryEnd, mTargetStart, mTargetEnd;
	sequence mAlignedQuery, mAlignedTarget;
	double mBitScore;
	double mExpect;
	bool mGapped;

	bool operator>(const BlastHsp &inHsp) const { return mScore > inHsp.mScore; }
	void CalculateExpect(int64_t inSearchSpace, double inLambda, double inLogKappa);

	bool Overlaps(const BlastHsp &inOther) const
	{
		return mQueryEnd >= inOther.mQueryStart and mQueryStart <= inOther.mQueryEnd and
		       mTargetEnd >= inOther.mTargetStart and mTargetStart <= inOther.mTargetEnd;
	}

	size_t length() const
	{
		return mAlignedQuery.length();
	}

	float identity()
	{
		assert(mAlignedQuery.length() == mAlignedTarget.length());

		size_t n = 0, L = mAlignedQuery.length();

		for (size_t i = 0; i < L; ++i)
		{
			if (mAlignedQuery[i] == mAlignedTarget[i])
				++n;
		}

		return static_cast<float>(n) / L;
	}
};

struct BlastHit
{
	std::string mDefLine;
	sequence mTarget;
	std::vector<BlastHsp> mHsps;

	BlastHit(const std::string &inDefLine, const sequence &inTarget)
		: mDefLine(inDefLine)
		, mTarget(inTarget)
	{
	}

	BlastHit(const BlastHit &hit) = default;
	BlastHit(BlastHit &&hit) = default;
};

// --------------------------------------------------------------------

/// \brief Search in \a inDatabank for \a inQuery using the blastp algorithm
///
/// \param inDatabank		The databank to search, should be a FastA file
/// \param inQuery			The protein sequence to use as query
/// \param inMatrix			The scoring matrix to use
/// \param inWordSize		The word size, values can be 2, 3 or 4
/// \param inExpect			The maximum expect value for hits, usually 10
/// \param inFilter			When true a low complexity filter is used
/// \param inGapped			When true, gapped alignment is done
/// \param inGapOpen		The cost of opening a gap
/// \param inGapExtend		The cost of extending a gap
/// \param inReportLimit	The maximum number of hits to return
/// \param inThreads		Maximum number of threads to use
/// \result					An array of blasthit structs

std::vector<BlastHit> BlastP(const std::filesystem::path &inDatabank,
	const std::string &inQuery, const std::string &inMatrix,
	uint32_t inWordSize, double inExpect, bool inFilter,
	bool inGapped, int32_t inGapOpen, int32_t inGapExtend,
	uint32_t inReportLimit, uint32_t inThreads);

/// \brief Search \a inDatabank for \a inQuery using the blastp algorithm, with default parameters
///
/// \param inDatabank	The databank to search, should be a FastA file
/// \param inQuery		The protein sequence to use as query
/// \result				An array of blasthit structs

inline std::vector<BlastHit> BlastP(const std::filesystem::path &inDatabank, const std::string &inQuery, uint32_t inReportLimit)
{
	return BlastP(inDatabank, inQuery, "BLOSUM62", 3, 10, true, true, 11, 1, inReportLimit, std::thread::hardware_concurrency());
}
