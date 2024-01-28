/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/ray_line_intersector.hpp"
#include "detray/navigation/intersection/ray_plane_intersector.hpp"

namespace detray {

/// @brief Intersection implementation for detector surfaces using a ray
/// trajectory.
///
/// @note specialized into the concrete intersectors for the different local
/// geometries in the respective header files
template <typename algebra_t, typename frame_t>
struct ray_intersector_impl {};

template <typename algebra_t, typename shape_t>
using ray_intersector = ray_intersector_impl<
    algebra_t, typename shape_t::template local_frame_type<algebra_t>>;

}  // namespace detray
