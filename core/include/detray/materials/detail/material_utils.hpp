/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/geometry/surface.hpp"

namespace detray::detail {

/// Helper function to calculate the incidence angle of a track on a surface
///
/// @param ctx the geometric context
/// @param sf the geometric surface
/// @param track_dir normalized direction vector of the track
/// @param loc the local/bound position of the track on the surface
///
/// @returns the cosine of the incidence angle given a local/bound position
template <typename detector_t, typename point_t>
    requires std::is_same_v<point_t,
                            dpoint3D<typename detector_t::algebra_type>> ||
    std::is_same_v<point_t, dpoint2D<typename detector_t::algebra_type>>
        DETRAY_HOST_DEVICE constexpr dscalar<typename detector_t::algebra_type>
        cos_angle(const typename detector_t::geometry_context &ctx,
                  geometry::surface<detector_t> sf,
                  const dvector3D<typename detector_t::algebra_type> &track_dir,
                  const point_t &loc) {
    return math::fabs(vector::dot(track_dir, sf.normal(ctx, loc)));
}

}  // namespace detray::detail
