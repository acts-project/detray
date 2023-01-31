/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracer/definitions/colors.hpp"

namespace detray::texture::detail {

/// @brief holds rgb and alpha values for color shading
template <typename depth = uint8_t>
DETRAY_HOST_DEVICE inline constexpr detray::texture::color<depth>
material_color_helper(const material<scalar> &mat) {
    // color based on material
    if (mat == detray::beryllium<scalar>{} or
        mat == detray::beryllium_tml<scalar>{}) {
        return detray::texture::grey<depth>;
    } else if (mat == detray::aluminium<scalar>{}) {
        return detray::texture::light_grey<depth>;
    } else if (mat == detray::tungsten<scalar>{}) {
        return detray::texture::dim_grey<depth>;
    } else if (mat == detray::gold<scalar>{}) {
        return detray::texture::golden_rod<depth>;
    } else if (mat == detray::silicon<scalar>{} or
               mat == detray::silicon_tml<scalar>{}) {
        return detray::texture::dark_grey<depth>;
    } else {
        // default for unknown material
        return detray::texture::dark_red<depth>;
    }
}

}  // namespace detray::texture::detail
