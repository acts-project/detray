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

namespace detray {

// Rod structure to be mapped on the mask
template <typename scalar_t>
struct material_rod {
    using scalar_type = scalar_t;
    using material_type = material<scalar_t>;

    material_rod() = default;

    material_rod(const material_type& material, scalar_type radius)
        : m_material(material), m_radius(radius) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    bool operator==(const material_rod<scalar_t>& rhs) const {
        return (m_material == rhs.get_material() && m_radius == rhs.radius());
    }

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_type& get_material() const { return m_material; }
    /// Return the radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type radius() const { return m_radius; }
    /// Return the pre interaction length
    DETRAY_HOST_DEVICE
    constexpr scalar_type pre_interaction_length(
        const line_plane_intersection& is) const {
        return m_radius * is.cos_incidence_angle;
    }
    /// Return the post interaction length
    DETRAY_HOST_DEVICE
    constexpr scalar_type post_interaction_length(
        const line_plane_intersection& is) const {
        return m_radius * is.cos_incidence_angle;
    }

    private:
    material_type m_material = {};
    scalar_type m_radius = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray