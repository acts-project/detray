/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/definitions/qualifiers.hpp"

// object registry for index geometry
template <typename surface_type, typename portal_type>
struct object_registry {
    // Known primitives
    enum id : unsigned int {
        e_object_types = 2,
        e_surface = 0,
        e_portal = 1,
        e_any = 1,  // defaults to portal
        e_unknown = 3,
    };

    template <typename value_type>
    DETRAY_HOST_DEVICE static constexpr auto get() {
        if constexpr (std::is_same_v<value_type, surface_type>) {
            return e_surface;
        }
        if constexpr (std::is_same_v<value_type, portal_type>) {
            return e_portal;
        }
        return e_unknown;
    }
};

// object registry for unified index geometry
template <typename surface_type>
struct unified_object_registry {
    // Known primitives
    enum id : unsigned int {
        e_object_types = 1,
        e_surface = 0,
        e_portal = 0,  // not used (same as surface)
        e_any = 0,
        e_unknown = 2,
    };

    template <typename value_type>
    DETRAY_HOST_DEVICE static constexpr auto get() {
        if constexpr (std::is_same_v<value_type, surface_type>) {
            return e_surface;
        }
        return e_unknown;
    }
};

// Types for toy geometry
struct toy_object_registry {
    // Known primitives
    enum id : unsigned int {
        e_object_types = 1,
        e_surface = 0,
        e_portal = 0,  // same as surface
        e_any = 0,
        e_unknown = 2,
    };

    template <typename value_type = void>
    DETRAY_HOST_DEVICE static constexpr auto get() {
        return e_surface;
    }
};
