/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/materials/detail/concepts.hpp"
#include "detray/materials/detail/material_accessor.hpp"
#include "detray/materials/material.hpp"

namespace detray::detail {

/// A functor to retrieve the volume material parameters at a given position
struct get_material_params {
    template <typename mat_group_t, typename index_t, concepts::point point_t>
    DETRAY_HOST_DEVICE inline auto operator()(const mat_group_t &mat_group,
                                              const index_t &idx,
                                              const point_t &loc_p) const {
        using material_t = typename mat_group_t::value_type;

        if constexpr (concepts::volume_material<material_t>) {

            if constexpr (concepts::homogeneous_material<material_t>) {
                // Homogeneous volume material
                return &(detail::material_accessor::get(mat_group, idx, loc_p));
            } else {
                // Volume material maps
                return &(detail::material_accessor::get(mat_group, idx, loc_p)
                             .get_material());
            }
        } else {
            using scalar_t = typename material_t::scalar_type;
            // Cannot be reached for volumes
            return static_cast<const material<scalar_t> *>(nullptr);
        }
    }
};

/// A functor to access all surfaces registered in an acceleration structure
/// @tparam functor_t functor that performs a task per surface
template <typename functor_t>
struct apply_to_surfaces {

    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename accel_group_t, typename accel_index_t, typename... Args>
    DETRAY_HOST_DEVICE inline void operator()(const accel_group_t &group,
                                              const accel_index_t index,
                                              Args &&...args) const {

        // Iterate through the surface neighbouhood of the track
        for (const auto &sf : group[index].all()) {
            functor_t{}(sf, std::forward<Args>(args)...);
        }
    }
};

/// A functor to find surfaces in the neighborhood of a track position
/// @tparam functor_t functor that performs a task per surface in the
///                   neighbourhood (e.g. intersection)
template <typename functor_t>
struct apply_to_neighbourhood {

    /// Call operator that forwards the neighborhood search call in a volume
    /// to a surface finder data structure
    template <typename accel_group_t, typename accel_index_t,
              typename detector_t, typename track_t, typename config_t,
              typename... Args>
    DETRAY_HOST_DEVICE inline void operator()(
        const accel_group_t &group, const accel_index_t index,
        const detector_t &det, const typename detector_t::volume_type &volume,
        const track_t &track, const config_t &cfg,
        const typename detector_t::geometry_context &ctx,
        Args &&...args) const {

        // Get the acceleration structure for the volume
        decltype(auto) accel = group[index];

        // Iterate through the surface neighbouhood of the track
        for (const auto &sf : accel.search(det, volume, track, cfg, ctx)) {
            functor_t{}(sf, std::forward<Args>(args)...);
        }
    }
};

}  // namespace detray::detail
