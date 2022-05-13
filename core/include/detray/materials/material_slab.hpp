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

template <typename material_t>
struct material_slab {
    using material_type = material_t;
    using scalar_type = typename material_t::scalar_type;

    material_slab(const material_t& material, scalar_type thickness)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()) {}

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_t& material() const { return m_material; }
    /// Return the thickness.
    DETRAY_HOST_DEVICE
    constexpr scalar_type thickness() const { return m_thickness; }
    /// Return the radiation length fraction.
    DETRAY_HOST_DEVICE
    constexpr scalar_type thickness_in_X0() const { return m_thickness_in_X0; }
    /// Return the nuclear interaction length fraction.
    DETRAY_HOST_DEVICE
    constexpr scalar_type thickness_in_L0() const { return m_thickness_in_L0; }

    private:
    material_t m_material;
    scalar_type m_thickness;
    scalar_type m_thickness_in_X0;
    scalar_type m_thickness_in_L0;
};

}  // namespace detray