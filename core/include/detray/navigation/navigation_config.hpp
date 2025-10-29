/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"

// System include(s)
#include <ostream>

namespace detray::navigation {

/// Navigation configuration
struct config {
    /// Tolerance on the mask 'is_inside' check:
    intersection::config intersection{};
    /// Search window size for grid based acceleration structures
    /// (0, 0): only look at current bin
    darray<dindex, 2> search_window = {0u, 0u};

    // Actor configuration

    /// Percentage of total track path to assume as accumulated error
    float accumulated_error{0.001f};
    /// Number of standard deviations to assume to model the scattering noise
    int n_scattering_stddev{2};
    /// Add adaptive mask tolerance to navigation
    bool estimate_scattering_noise{true};

    /// Print the navigation configuration
    DETRAY_HOST
    friend std::ostream& operator<<(std::ostream& out, const config& cfg) {
        out << cfg.intersection
            << "  Search window         : " << cfg.search_window[0] << " x "
            << cfg.search_window[1] << "\n";

        if (cfg.estimate_scattering_noise) {
            out << "Actor configuration:\n"
                << "  Accumulated error     : " << cfg.accumulated_error * 100.f
                << " %\n"
                << "  No. scattering stddev : " << cfg.n_scattering_stddev
                << "\n";
        }

        return out;
    }
};
}  // namespace detray::navigation
