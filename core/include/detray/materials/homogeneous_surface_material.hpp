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

template <typename material_structure_t>
struct homogeneous_surface_material {

    using scalar_type = typename material_structure_t::scalar_type;
    using point2 = __plugin::point2<scalar_t>;
    using point3 = __plugin::point3<scalar_t>;

    homogeneous_surface_material(const material_structure_t& structure)
        : m_structure(structure) {}

    DETRAY_HOST_DEVICE
    const material_structure_t& material_structure(const point2& /*p*/) const {
        return m_structure;
    }

    DETRAY_HOST_DEVICE
    const material_structure_t& material_structure(const point3& /*p*/) const {
        return m_structure;
    }

    private:
    material_structure_t m_structure;
};

}  // namespace detray