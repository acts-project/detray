/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::detail {

/// A functor to access the surfaces of a volume
template <typename functor_t>
struct surface_getter {

    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename accel_group_t, typename accel_index_t, typename... Args>
    DETRAY_HOST_DEVICE inline void operator()(const accel_group_t &group,
                                              const accel_index_t index,
                                              Args &&... args) const {

        // Run over the surfaces in a single acceleration data structure
        for (const auto &sf : group[index].all()) {
            functor_t{}(sf, std::forward<Args>(args)...);
        }
    }
};

/// A functor to find surfaces in the neighborhood of a track position
template <typename functor_t>
struct neighborhood_getter {

    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename accel_group_t, typename accel_index_t,
              typename detector_t, typename track_t, typename config_t,
              typename... Args>
    DETRAY_HOST_DEVICE inline void operator()(
        const accel_group_t &group, const accel_index_t index,
        const detector_t &det, const typename detector_t::volume_type &volume,
        const track_t &track, const config_t &cfg, Args &&... args) const {

        decltype(auto) accel = group[index];

        // Run over the surfaces in a single acceleration data structure
        for (const auto &sf : accel.search(det, volume, track, cfg)) {
            functor_t{}(sf, std::forward<Args>(args)...);
        }
    }
};

}  // namespace detray::detail
