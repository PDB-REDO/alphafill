// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at
//             http://www.boost.org/LICENSE_1_0.txt)
//
// substitution matrix for multiple sequence alignments

#pragma once

// #include "mas.h"

#include <algorithm>
#include <cassert>
#include <istream>
#include <stdexcept>
#include <string>
#include <vector>

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

// --------------------------------------------------------------------
// uBlas compatible matrix types
// matrix is m x n, addressing i,j is 0 <= i < m and 0 <= j < n
// element m i,j is mapped to [i * n + j] and thus storage is row major

template <typename T>
class matrix_base
{
  public:
	typedef T value_type;

	virtual ~matrix_base() {}

	virtual uint32_t dim_m() const = 0;
	virtual uint32_t dim_n() const = 0;

	virtual value_type &operator()(uint32_t i, uint32_t j) { throw std::runtime_error("unimplemented method"); }
	virtual value_type operator()(uint32_t i, uint32_t j) const = 0;

	matrix_base &operator*=(const value_type &rhs);

	matrix_base &operator-=(const value_type &rhs);
};

template <typename T>
matrix_base<T> &matrix_base<T>::operator*=(const T &rhs)
{
	for (uint32_t i = 0; i < dim_m(); ++i)
	{
		for (uint32_t j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) *= rhs;
		}
	}

	return *this;
}

template <typename T>
matrix_base<T> &matrix_base<T>::operator-=(const T &rhs)
{
	for (uint32_t i = 0; i < dim_m(); ++i)
	{
		for (uint32_t j = 0; j < dim_n(); ++j)
		{
			operator()(i, j) -= rhs;
		}
	}

	return *this;
}

template <typename T>
std::ostream &operator<<(std::ostream &lhs, const matrix_base<T> &rhs)
{
	lhs << '[' << rhs.dim_m() << ',' << rhs.dim_n() << ']' << '(';
	for (uint32_t i = 0; i < rhs.dim_m(); ++i)
	{
		lhs << '(';
		for (uint32_t j = 0; j < rhs.dim_n(); ++j)
		{
			if (j > 0)
				lhs << ',';
			lhs << rhs(i, j);
		}
		lhs << ')';
	}
	lhs << ')';

	return lhs;
}

template <typename T>
class matrix : public matrix_base<T>
{
  public:
	typedef T value_type;

	template <typename T2>
	matrix(const matrix_base<T2> &m)
		: m_m(m.dim_m())
		, m_n(m.dim_n())
	{
		m_data = new value_type[m_m * m_n];
		for (uint32_t i = 0; i < m_m; ++i)
		{
			for (uint32_t j = 0; j < m_n; ++j)
				operator()(i, j) = m(i, j);
		}
	}

	matrix()
		: m_data(nullptr)
		, m_m(0)
		, m_n(0)
	{
	}

	matrix(const matrix &m)
		: m_m(m.m_m)
		, m_n(m.m_n)
	{
		m_data = new value_type[m_m * m_n];
		std::copy(m.m_data, m.m_data + (m_m * m_n), m_data);
	}

	matrix &operator=(const matrix &m)
	{
		value_type *t = new value_type[m.m_m * m.m_n];
		std::copy(m.m_data, m.m_data + (m.m_m * m.m_n), t);

		delete[] m_data;
		m_data = t;
		m_m = m.m_m;
		m_n = m.m_n;

		return *this;
	}

	matrix(uint32_t m, uint32_t n, T v = T())
		: m_m(m)
		, m_n(n)
	{
		m_data = new value_type[m_m * m_n];
		std::fill(m_data, m_data + (m_m * m_n), v);
	}

	virtual ~matrix()
	{
		delete[] m_data;
	}

	virtual uint32_t dim_m() const { return m_m; }
	virtual uint32_t dim_n() const { return m_n; }

	virtual value_type operator()(uint32_t i, uint32_t j) const
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	virtual value_type &operator()(uint32_t i, uint32_t j)
	{
		assert(i < m_m);
		assert(j < m_n);
		return m_data[i * m_n + j];
	}

	template <typename Func>
	void each(Func f)
	{
		for (uint32_t i = 0; i < m_m * m_n; ++i)
			f(m_data[i]);
	}

	template <typename U>
	matrix &operator/=(U v)
	{
		for (uint32_t i = 0; i < m_m * m_n; ++i)
			m_data[i] /= v;

		return *this;
	}

  private:
	value_type *m_data;
	uint32_t m_m, m_n;
};

// --------------------------------------------------------------------

template <typename T>
class symmetric_matrix : public matrix_base<T>
{
  public:
	typedef typename matrix_base<T>::value_type value_type;

	symmetric_matrix(uint32_t n, T v = T())
		: m_owner(true)
		, m_n(n)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;
		m_data = new value_type[N];
		std::fill(m_data, m_data + N, v);
	}

	symmetric_matrix(const T *data, uint32_t n)
		: m_owner(false)
		, m_data(const_cast<T *>(data))
		, m_n(n)
	{
	}

	virtual ~symmetric_matrix()
	{
		if (m_owner)
			delete[] m_data;
	}

	virtual uint32_t dim_m() const { return m_n; }
	virtual uint32_t dim_n() const { return m_n; }

	T operator()(uint32_t i, uint32_t j) const;
	virtual T &operator()(uint32_t i, uint32_t j);

	// erase two rows, add one at the end (for neighbour joining)
	void erase_2(uint32_t i, uint32_t j);

	template <typename Func>
	void each(Func f)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;

		for (uint32_t i = 0; i < N; ++i)
			f(m_data[i]);
	}

	template <typename U>
	symmetric_matrix &operator/=(U v)
	{
		uint32_t N = (m_n * (m_n + 1)) / 2;

		for (uint32_t i = 0; i < N; ++i)
			m_data[i] /= v;

		return *this;
	}

  private:
	bool m_owner;
	value_type *m_data;
	uint32_t m_n;
};

template <typename T>
inline T symmetric_matrix<T>::operator()(uint32_t i, uint32_t j) const
{
	return i < j
	           ? m_data[(j * (j + 1)) / 2 + i]
	           : m_data[(i * (i + 1)) / 2 + j];
	//	if (i > j)
	//		std::swap(i, j);
	//	assert(j < m_n);
	//	return m_data[(j * (j + 1)) / 2 + i];
}

template <typename T>
inline T &symmetric_matrix<T>::operator()(uint32_t i, uint32_t j)
{
	if (i > j)
		std::swap(i, j);
	assert(j < m_n);
	return m_data[(j * (j + 1)) / 2 + i];
}

template <typename T>
void symmetric_matrix<T>::erase_2(uint32_t di, uint32_t dj)
{
	uint32_t s = 0, d = 0;
	for (uint32_t i = 0; i < m_n; ++i)
	{
		for (uint32_t j = 0; j < i; ++j)
		{
			if (i != di and j != dj and i != dj and j != di)
			{
				if (s != d)
					m_data[d] = m_data[s];
				++d;
			}

			++s;
		}
	}

	--m_n;
}

template <typename T>
class identity_matrix : public matrix_base<T>
{
  public:
	typedef typename matrix_base<T>::value_type value_type;

	identity_matrix(uint32_t n)
		: m_n(n)
	{
	}

	virtual uint32_t dim_m() const { return m_n; }
	virtual uint32_t dim_n() const { return m_n; }

	virtual value_type operator()(uint32_t i, uint32_t j) const
	{
		value_type result = 0;
		if (i == j)
			result = 1;
		return result;
	}

  private:
	uint32_t m_n;
};

// --------------------------------------------------------------------
// matrix functions

template <typename T>
matrix<T> operator*(const matrix_base<T> &lhs, const matrix_base<T> &rhs)
{
	matrix<T> result(std::min(lhs.dim_m(), rhs.dim_m()), std::min(lhs.dim_n(), rhs.dim_n()));

	for (uint32_t i = 0; i < result.dim_m(); ++i)
	{
		for (uint32_t j = 0; j < result.dim_n(); ++j)
		{
			for (uint32_t li = 0, rj = 0; li < lhs.dim_m() and rj < rhs.dim_n(); ++li, ++rj)
				result(i, j) += lhs(li, j) * rhs(i, rj);
		}
	}

	return result;
}

template <typename T>
matrix<T> operator*(const matrix_base<T> &lhs, T rhs)
{
	matrix<T> result(lhs);
	result *= rhs;

	return result;
}

template <typename T>
matrix<T> operator-(const matrix_base<T> &lhs, const matrix_base<T> &rhs)
{
	matrix<T> result(std::min(lhs.dim_m(), rhs.dim_m()), std::min(lhs.dim_n(), rhs.dim_n()));

	for (uint32_t i = 0; i < result.dim_m(); ++i)
	{
		for (uint32_t j = 0; j < result.dim_n(); ++j)
		{
			result(i, j) = lhs(i, j) - rhs(i, j);
		}
	}

	return result;
}

template <typename T>
matrix<T> operator-(const matrix_base<T> &lhs, T rhs)
{
	matrix<T> result(lhs.dim_m(), lhs.dim_n());
	result -= rhs;
	return result;
}

template <typename T>
symmetric_matrix<T> hamming_distance(const matrix_base<T> &lhs, T rhs);

template <typename T>
std::vector<T> sum(const matrix_base<T> &m);

// #include "matrix.inl"
