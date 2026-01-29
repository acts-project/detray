/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/utils/concepts.hpp"

namespace detray {

/// Size around a query point for which to do the search in an acceleration
/// structure
///
/// @tparam window_size_t type that defines the window size (e.g. distance or
/// number of bins)
/// @tparam DIM number of dimensions of the window boundaries
template <concepts::arithmetic window_size_t, std::size_t DIM>
using search_window = darray<window_size_t, DIM>;

}  // namespace detray
