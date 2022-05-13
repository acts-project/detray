/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Detray include(s)
#include "detray/definitions/qualifiers.hpp"

namespace detray {

template <typename material_shape_t>
struct homogeneous_surface_material {

    using scalar_type = typename material_shape_t::scalar_type;
    using point2 = __plugin::point2<scalar_t>;
    using point3 = __plugin::point3<scalar_t>;

    homogeneous_surface_material(const material_shape_t& shape)
        : m_shape(shape) {}

    DETRAY_HOST_DEVICE
    const material_shape_t& material_shape(const point2& /*p*/) const {
        return m_shape;
    }

    DETRAY_HOST_DEVICE
    const material_shape_t& material_shape(const point3& /*p*/) const {
        return m_shape;
    }

    private:
    material_shape_t m_shape;
};

}  // namespace detray