// SPDX-PackageName: "detray, a part of the ACTS project"
// SPDX-FileCopyrightText: 2021 CERN
// SPDX-License-Identifier: MPL-2.0

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/intersection/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/ray_line_intersector.hpp"
#include "detray/navigation/intersection/ray_plane_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_cylinder_portal_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_line_intersector.hpp"
#include "detray/navigation/intersection/soa/ray_plane_intersector.hpp"

namespace detray {

/// @brief Intersection implementation for detector surfaces using a ray
/// trajectory.
///
/// @note specialized into the concrete intersectors for the different local
/// geometries in the respective header files
template <typename frame_t, concepts::algebra algebra_t, bool do_debug>
struct ray_intersector_impl {};

template <typename shape_t, concepts::algebra algebra_t, bool do_debug = false>
using ray_intersector =
    ray_intersector_impl<typename shape_t::template local_frame_type<algebra_t>,
                         algebra_t, do_debug>;

}  // namespace detray
