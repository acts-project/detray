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
    using transform3 = __plugin::transform3<scalar_t>;

    material_slab() = default;

    /// Constructor for material slab
    /// @param material is the elemental or mixture material
    /// @param thickness is the thickness of slab
    material_slab(const material_type& material, const scalar_type thickness)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()) {}

    /// Constructor for material slab
    /// @param material is the elemental or mixture material
    /// @param thickness is the thickness of slab
    /// @param interaction_point is the interaction point where the measurement
    /// happens
    material_slab(const material_type& material, const scalar_type thickness,
                  const scalar_type interaction_point)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()),
          m_interaction_point(interaction_point) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE bool operator==(
        const material_slab<scalar_t>& rhs) const {
        return (m_material == rhs.get_material() &&
                m_thickness == rhs.thickness());
    }

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

    // interaction point is the point where the measurement happens in the slab
    // 0 (epsilon) means the middle of the slab
    scalar_type m_interaction_point = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray