/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/detector_volume.hpp"

namespace detray::svgtools::utils {
/// @return the sub-surfaces of the volume - non-const access
template <typename detector_t>
auto surface_lookup(const detector_t& detector,
                    const detector_volume<detector_t>& d_volume) {
    typename detector_t::surface_lookup_container descriptors;
    for (auto desc : detector.surface_lookup()) {
        if (desc.volume() == d_volume.index()) {
            descriptors.push_back(desc);
        }
    }
    return descriptors;
}
}  // namespace detray::svgtools::utils
