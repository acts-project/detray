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

// Rod structure to be mapped on the mask
template <typename material_t>
struct material_rod {
    using material_type = material_t;
    using scalar_type = typename material_t::scalar_type;

    material_rod(const material_t& material, scalar_type outer_r,
                 scalar_type inner_r = 0)
        : m_material(material), m_outer_r(outer_r), m_inner_r(inner_r) {}

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_t& material() const { return m_material; }
    /// Return the outer radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type outer_r() const { return m_outer_r; }
    /// Return the inner radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type inner_r() const { return m_inner_r; }

    private:
    material_t m_material;
    scalar_type m_outer_r;
    scalar_type m_inner_r;
};

}  // namespace detray