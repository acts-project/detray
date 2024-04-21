/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

// System include(s)
#include <ostream>

namespace detray::navigation {

/// Navigation trust levels determine how the candidates chache is updated
enum class trust_level {
    e_no_trust = 0,  ///< re-initialize the volume (i.e. run local navigation)
    e_fair = 1,      ///< update the distance & order of the candidates
    e_high = 3,  ///< update the distance to the next candidate (current target)
    e_full = 4   ///< don't update anything
};

/// Navigation configuration
struct config {
    /// Tolerance on the mask 'is_inside' check:
    /// @{
    /// Minimal tolerance: ~ position uncertainty on surface
    float min_mask_tolerance{1e-5f * unit<float>::mm};
    /// Maximal tolerance: loose tolerance when still far away from surface
    float max_mask_tolerance{1.f * unit<float>::mm};
    ///@}
    /// Maximal absolute path distance for a track to be considered 'on surface'
    float path_tolerance{1.f * unit<float>::um};
    /// How far behind the track position to look for candidates
    float overstep_tolerance{-100.f * unit<float>::um};
    /// Search window size for grid based acceleration structures
    /// (0, 0): only look at current bin
    std::array<dindex, 2> search_window = {0u, 0u};
};

/// Print the navigation configuration
DETRAY_HOST
inline std::ostream& operator<<(std::ostream& out,
                                const detray::navigation::config& cfg) {
    out << "  Minimal mask tolerance: "
        << cfg.min_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
        << "  Maximal mask tolerance: "
        << cfg.max_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
        << "  Path tolerance        : "
        << cfg.path_tolerance / detray::unit<float>::um << " [um]\n"
        << "  Overstep tolerance    : "
        << cfg.overstep_tolerance / detray::unit<float>::um << " [um]\n"
        << "  Search window         : " << cfg.search_window[0] << " x "
        << cfg.search_window[1] << "\n";

    return out;
}

}  // namespace detray::navigation
