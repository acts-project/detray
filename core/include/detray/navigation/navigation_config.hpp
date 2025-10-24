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

// System include(s)
#include <ostream>

namespace detray::navigation {

/// Navigation configuration
struct config {
    /// Tolerance on the mask 'is_inside' check:
    /// @{
    /// Minimal tolerance: ~ position uncertainty on surface
    float min_mask_tolerance{1e-5f * unit<float>::mm};
    /// Maximal tolerance: loose tolerance when still far away from surface
    float max_mask_tolerance{3.f * unit<float>::mm};
    /// Scale factor on the path used for the mask tolerance calculation
    float mask_tolerance_scalor{5e-2f};
    /// @}
    /// Maximal absolute path distance for a track to be considered 'on surface'
    float path_tolerance{1.f * unit<float>::um};
    /// How far behind the track position to look for candidates
    float overstep_tolerance{-1000.f * unit<float>::um};
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
        out << "  Min. mask tolerance   : "
            << cfg.min_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
            << "  Max. mask tolerance   : "
            << cfg.max_mask_tolerance / detray::unit<float>::mm << " [mm]\n"
            << "  Mask tolerance scalor : " << cfg.mask_tolerance_scalor << "\n"
            << "  Path tolerance        : "
            << cfg.path_tolerance / detray::unit<float>::um << " [um]\n"
            << "  Overstep tolerance    : "
            << cfg.overstep_tolerance / detray::unit<float>::um << " [um]\n"
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
