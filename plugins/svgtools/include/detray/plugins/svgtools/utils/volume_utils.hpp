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

template <typename detector_t, typename dindex>
auto surface_indices(const detector_t& detector,
                    const dindex volume_index) {
     const auto d_volume =
            detector.volume_by_index(static_cast<detray::dindex>(volume_index));
    const auto descriptors = surface_lookup(detector, d_volume);
    std::vector<dindex> ret;
    for (const auto& desc : descriptors){
        ret.push_back(desc.index());
    }
    return ret;
}

}  // namespace detray::svgtools::utils
