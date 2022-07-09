/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
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

    /// Constructor
    /// @param material is the elemental or mixture material
    /// @param thickness is the thickness of the slab
    material_slab(const material_type& material, scalar_type thickness)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()) {}

    /// Constructor
    /// @param material is the elemental or mixture material
    /// @param thickness is the thickness of the slab
    /// @param interaction_point is the interaction point of the line
    material_slab(const material_type& material, scalar_type thickness,
                  scalar_type interaction_point)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()),
          m_interaction_point(interaction_point) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    bool operator==(const material_slab<scalar_t>& rhs) const {
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
    /// Return the pre interaction length
    DETRAY_HOST_DEVICE
    constexpr scalar_type pre_interaction_length(
        const line_plane_intersection& is) const {
        if (is.direction == intersection::direction::e_along) {
            return (m_thickness / scalar_type(2) + m_interaction_point) /
                   is.cos_incidence_angle;
        } else if (is.direction == intersection::direction::e_opposite) {
            return (m_thickness / scalar_type(2) - m_interaction_point) /
                   is.cos_incidence_angle;
        } else {
            return 0;
        }
    }
    /// Return the post interaction length
    DETRAY_HOST_DEVICE
    constexpr scalar_type post_interaction_length(
        const line_plane_intersection& is) const {
        if (is.direction == intersection::direction::e_along) {
            return (m_thickness / scalar_type(2) - m_interaction_point) /
                   is.cos_incidence_angle;
        } else if (is.direction == intersection::direction::e_opposite) {
            return (m_thickness / scalar_type(2) + m_interaction_point) /
                   is.cos_incidence_angle;
        } else {
            return 0;
        }
    }

    private:
    material_type m_material = {};
    scalar_type m_thickness = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_L0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_interaction_point = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray