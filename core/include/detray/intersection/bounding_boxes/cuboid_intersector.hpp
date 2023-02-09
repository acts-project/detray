/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"

// System include(s)
#include <cmath>
#include <limits>

namespace detray {

/// Intersector that checks if a cuboid bonuding box is intersected by a ray.
///
/// This will not compute the intersection point and is meant for aabb
/// intersection instead of navigation.
/// For the implementation see:
/// https://tavianator.com/2022/ray_box_boundary.html
/// and
/// https://github.com/acts-project/acts/blob/85a67758f79f3b01aed4025592ab6a5e5ffbd323/Core/include/Acts/Utilities/BoundingBox.ipp
struct cuboid_intersector {

    /// Operator to find intersections between ray and a cuboid aabb
    ///
    /// @tparam algebra_t is the linear algebra implementation used
    ///
    /// @param ray is the input ray trajectory
    /// @param box is the input mask
    /// @param mask_tolerance is the tolerance for mask edges
    ///
    /// @return true if the cuboid aabb was intersected
    template <typename algebra_t, typename mask_t>
    DETRAY_HOST_DEVICE bool operator()(
        const detail::ray<algebra_t> &ray, const mask_t &box,
        const typename algebra_t::scalar_type /*mask_tolerance*/ = 0.f) const {

        // using scalar_t = typename algebra_t::scalar_type;
        using point3 = typename algebra_t::point3;
        using vector3 = typename algebra_t::vector3;
        using boundaries = typename mask_t::boundaries;

        const point3 &ro = ray.pos();
        const vector3 &rd = ray.dir();
        // @todo put operator/ in algebra-plugins
        const vector3 inv_dir{1.f / rd[0], 1.f / rd[1], 1.f / rd[2]};

        // This is prob. slow -> @todo refactor masks to hold custom mask values
        const point3 min{box[boundaries::e_min_x], box[boundaries::e_min_y],
                         box[boundaries::e_min_z]};
        const point3 max{box[boundaries::e_max_x], box[boundaries::e_max_y],
                         box[boundaries::e_max_z]};
        // const auto* min = new (box.values().data()) point3();
        // const auto* max = new (box.values().data() + 3) point3();

        // scalar_t tmin{0.f}, tmax{std::numeric_limits<scalar_t>::infinity()};

        const point3 tmin = (min - ro) * inv_dir;
        const point3 tmax = (max - ro) * inv_dir;

        return true;
    }
};

}  // namespace detray