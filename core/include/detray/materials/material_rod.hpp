/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/qualifiers.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"

namespace detray {

// Rod structure to be mapped on the line mask
template <typename scalar_t>
struct material_rod {
    using scalar_type = scalar_t;
    using material_type = material<scalar_t>;

    constexpr material_rod() = default;

    constexpr material_rod(const material_type& material, scalar_type radius)
        : m_material(material),
          m_radius(radius),
          m_radius_in_X0(radius / material.X0()),
          m_radius_in_L0(radius / material.L0()) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE
    constexpr bool operator==(const material_rod<scalar_t>& rhs) const {
        return (m_material == rhs.get_material() && m_radius == rhs.radius());
    }

    /// Boolean operator
    DETRAY_HOST_DEVICE
    constexpr operator bool() const {
        if (m_radius <= std::numeric_limits<scalar>::epsilon() ||
            m_material == vacuum<scalar_type>() || m_material.Z() == 0 ||
            m_material.mass_density() == 0) {
            return false;
        }
        return true;
    }

    /// Access the (average) material parameters.
    DETRAY_HOST_DEVICE
    constexpr const material_type& get_material() const { return m_material; }

    /// Return the radius
    DETRAY_HOST_DEVICE
    constexpr scalar_type radius() const { return m_radius; }

    /// Return the path segment
    DETRAY_HOST_DEVICE
    scalar_type path_segment(const line_plane_intersection& is) const {
        // Assume that is.p2[0] is radial distance of line intersector
        if (is.p2[0] > m_radius) {
            return 0.f;
        }

        const scalar_type sin_incidence_angle_2{
            1.f - is.cos_incidence_angle * is.cos_incidence_angle};

        return 2.f * std::sqrt((m_radius * m_radius - is.p2[0] * is.p2[0]) /
                               sin_incidence_angle_2);
    }
    /// Return the path segment in X0
    DETRAY_HOST_DEVICE
    scalar_type path_segment_in_X0(const line_plane_intersection& is) const {
        return this->path_segment(is) / m_material.X0();
    }
    /// Return the path segment in L0
    DETRAY_HOST_DEVICE
    scalar_type path_segment_in_L0(const line_plane_intersection& is) const {
        return this->path_segment(is) / m_material.L0();
    }

    private:
    material_type m_material = {};
    scalar_type m_radius = std::numeric_limits<scalar>::epsilon();
    scalar_type m_radius_in_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_radius_in_L0 = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray