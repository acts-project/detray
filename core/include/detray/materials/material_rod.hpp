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

namespace detray {

// Rod structure to be mapped on the mask
template <typename scalar_t>
struct material_rod {
    using scalar_type = scalar_t;
    using material_type = material<scalar_t>;

    material_rod() = default;

    material_rod(const material_type& material, scalar_type outer_r,
                 scalar_type inner_r = 0)
        : m_material(material), m_outer_r(outer_r), m_inner_r(inner_r) {}

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_type& get_material() const { return m_material; }
    /// Return the outer radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type outer_r() const { return m_outer_r; }
    /// Return the inner radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type inner_r() const { return m_inner_r; }

    private:
    material_type m_material = {};
    scalar_type m_outer_r = std::numeric_limits<scalar>::epsilon();
    scalar_type m_inner_r = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray