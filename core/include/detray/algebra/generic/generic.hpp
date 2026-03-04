/** Algebra plugins library, part of the ACTS project
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Impl include(s).
#include "detray/algebra/boolean.hpp"
#include "detray/algebra/impl/generic_matrix.hpp"
#include "detray/algebra/impl/generic_transform3.hpp"
#include "detray/algebra/impl/generic_vector.hpp"
#include "detray/algebra/math.hpp"

// Algorithms include(s).
#include "detray/algebra/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "detray/algebra/algorithms/matrix/determinant/cofactor.hpp"
#include "detray/algebra/algorithms/matrix/determinant/hard_coded.hpp"
#include "detray/algebra/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "detray/algebra/algorithms/matrix/inverse/cofactor.hpp"
#include "detray/algebra/algorithms/matrix/inverse/hard_coded.hpp"
#include "detray/algebra/algorithms/matrix/inverse/partial_pivot_lud.hpp"
#include "detray/algebra/algorithms/utils/algorithm_selector.hpp"
