/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Impl include(s).
#include "algebra/math/boolean.hpp"
#include "algebra/math/common.hpp"
#include "algebra/math/impl/generic_matrix.hpp"
#include "algebra/math/impl/generic_transform3.hpp"
#include "algebra/math/impl/generic_vector.hpp"

// Algorithms include(s).
#include "algebra/math/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "algebra/math/algorithms/matrix/determinant/cofactor.hpp"
#include "algebra/math/algorithms/matrix/determinant/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "algebra/math/algorithms/matrix/inverse/cofactor.hpp"
#include "algebra/math/algorithms/matrix/inverse/hard_coded.hpp"
#include "algebra/math/algorithms/matrix/inverse/partial_pivot_lud.hpp"
#include "algebra/math/algorithms/utils/algorithm_finder.hpp"
