/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/qualifiers.hpp"

// object registry for unified index geometry
template <typename surface_t>
struct object_registry {
    // Known primitives
    enum id : unsigned int {
        e_object_types = 1,
        e_surface = 0,
        e_portal = 0,  // not used (same as surface)
        e_any = 0,
        e_unknown = 2,
    };

    template <typename value_t>
    DETRAY_HOST_DEVICE static constexpr auto get() {
        if constexpr (std::is_same_v<value_t, surface_t>) {
            return e_surface;
        }
        return e_unknown;
    }
};