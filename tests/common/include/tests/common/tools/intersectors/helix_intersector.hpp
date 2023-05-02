/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/// @brief Intersection implementation for detector surfaces using a helix
/// trajectory.
///
/// The algorithm uses the Newton-Raphson method to find an intersection on
/// the unbounded surface and then applies the mask.
///
/// @note specialized into @c helix_plane_intersector and
/// @c helix_cylinder_intersector
template <typename intersection_t, typename mask_t, typename = void>
struct helix_intersector {};

}  // namespace detray