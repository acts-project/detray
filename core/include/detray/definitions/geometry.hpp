/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <cstdint>
#include <detray/definitions/detail/qualifiers.hpp>
#include <ostream>

namespace detray {

/// Shape of a detector volume.
///
/// cylinder: cylinder and two disc portals (shape: cylinder2D or cylinder3D).
/// cone: a cone portal and a disc portal (shape: missing).
/// rectangle: six rectangular portals that form a box (shape: cuboid3D).
/// trapezoid: six trapezoid portals (shape: cuboid3D).
/// cuboid: general cuboid form, excluding the previous ones (shape: cuboid3D).
enum class volume_id : std::uint_least8_t {
    e_cylinder = 0u,
    e_rectangle = 1u,
    e_trapezoid = 2u,
    e_cone = 3u,
    e_cuboid = 4u,
    e_unknown = 5u
};

#define _enum_print(x) \
    case x:            \
        os << #x;      \
        break

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, volume_id vid) {

    switch (vid) {
        using enum volume_id;
        _enum_print(e_cylinder);
        _enum_print(e_rectangle);
        _enum_print(e_trapezoid);
        _enum_print(e_cone);
        _enum_print(e_cuboid);
        _enum_print(e_unknown);
    }
    return os;
}

/// surface type, resolved during navigation.
///
/// sensitive: can provide measurements and have material.
/// passive: no measurements, but can have material.
/// portal: boundary surface between two detector volumes, can have material.
enum class surface_id : std::uint_least8_t {
    e_portal = 0u,
    e_sensitive = 1u,
    e_passive = 2u,
    e_unknown = 3u,
    e_all = e_unknown
};

DETRAY_HOST inline std::ostream& operator<<(std::ostream& os, surface_id sid) {

    switch (sid) {
        using enum surface_id;
        _enum_print(e_portal);
        _enum_print(e_sensitive);
        _enum_print(e_passive);
        _enum_print(e_unknown);
    }
    return os;
}

#undef _enum_print
}  // namespace detray
