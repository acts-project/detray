/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"

namespace detray {

/// A functor to find surfaces in the neighborhood of a track position
struct neighborhood_getter {

    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename sf_finder_group_t, typename sf_finder_index_t,
              typename detector_t, typename track_t>
    DETRAY_HOST_DEVICE inline auto operator()(
        const sf_finder_group_t &group, const sf_finder_index_t index,
        const detector_t &detector,
        const typename detector_t::volume_type &volume,
        const track_t &track) const {

        // Get surface finder for volume and perform the surface neighborhood
        // lookup
        const auto &sf_finder = group[index];
        return sf_finder.search(detector, volume, track);
    }
};

}  // namespace detray