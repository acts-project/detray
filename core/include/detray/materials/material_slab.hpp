/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/materials/material.hpp"

// System include(s)
#include <climits>


namespace detray {

// Slab structure to be mapped on the mask
template <typename scalar_t>
struct material_slab {
    using scalar_type = scalar_t;
    using material_type = material<scalar_t>;

    material_slab() = default;

    material_slab(const material_type& material, scalar_type thickness)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()) {}

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_type& get_material() const { return m_material; }
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
    material_type m_material = {};
    scalar_type m_thickness = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_L0 = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray