/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/units.hpp"

namespace detray::navigation {

/// Navigation trust levels determine how the candidates chache is updated
enum class trust_level {
    e_no_trust = 0,  ///< re-initialize the volume (i.e. run local navigation)
    e_fair = 1,      ///< update the distance & order of the candidates
    e_high = 3,  ///< update the distance to the next candidate (current target)
    e_full = 4   ///< don't update anything
};

/// Navigation configuration
template <typename scalar_t>
struct config {
    /// Tolerance on the masks 'is_inside' check
    scalar_t mask_tolerance{15.f * unit<scalar_t>::um};
    /// Maximal absolute path distance for a track to be considered 'on surface'
    scalar_t on_surface_tolerance{1.f * unit<scalar_t>::um};
    /// How far behind the track position to look for candidates
    scalar_t overstep_tolerance{-100.f * unit<scalar_t>::um};
    /// Search window size for grid based acceleration structures
    std::array<dindex, 2> search_window = {0u, 0u};
};

}  // namespace detray::navigation
