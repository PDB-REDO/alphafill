// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// align-3d

#pragma once

#include <cif++/Point.hpp>

class substitution_matrix_family;

double CalculateRMSD(const std::vector<mmcif::Point>& structureA, const std::vector<mmcif::Point>& structureB);

void AlignIterative(std::vector<mmcif::Point>& cAlphaA, std::vector<mmcif::Point>& cAlphaB,
	mmcif::Point& outTranslationA, mmcif::Point& outTranslationB, mmcif::quaternion& outRotation);
