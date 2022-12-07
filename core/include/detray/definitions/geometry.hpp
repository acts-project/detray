/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

/// Shape of a detector volume.
///
/// cylinder: cylinder and two disc portals (shape: cylinder2D or cylinder3D).
/// cone: a cone portal and a disc portal (shape: missing).
/// rectangle: six rectangular portals that form a box (shape: cuboid3D).
/// trapezoid: six trapezoid portals (shape: cuboid3D).
/// cuboid: general cuboid form, excluding the previous ones (shape: cuboid3D).
enum class volume_id {
    e_cylinder = 0,
    e_rectangle = 1,
    e_trapezoid = 2,
    e_cone = 3,
    e_cuboid = 4
};

/// surface type, resolved during navigation.
///
/// sensitive: can provide measurements and have material.
/// passive: no measurements, but can have material.
/// portal: boundary surface between two detector volumes, can have material.
enum class surface_id {
    e_sensitive = 0,
    e_portal = 1,
    e_passive = 2,
};

}  // namespace detray