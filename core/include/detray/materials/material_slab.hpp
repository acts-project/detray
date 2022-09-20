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
#include "detray/materials/predefined_materials.hpp"

// System include(s)
#include <climits>

namespace detray {

// Slab structure to be mapped on the mask (plane, cylinder)
template <typename scalar_t>
struct material_slab {
    using scalar_type = scalar_t;
    using material_type = material<scalar_t>;
    using transform3 = __plugin::transform3<scalar_t>;

    material_slab() = default;

    /// Constructor
    /// @param material is the elemental or mixture material
    /// @param thickness is the thickness of the slab
    material_slab(const material_type& material, scalar_type thickness)
        : m_material(material),
          m_thickness(thickness),
          m_thickness_in_X0(thickness / material.X0()),
          m_thickness_in_L0(thickness / material.L0()) {}

    /// Equality operator
    ///
    /// @param rhs is the right hand side to be compared to
    DETRAY_HOST_DEVICE bool operator==(
        const material_slab<scalar_t>& rhs) const {
        return (m_material == rhs.get_material() &&
                m_thickness == rhs.thickness());
    }

    /// Boolean operator
    DETRAY_HOST_DEVICE constexpr operator bool() const {
        if (m_thickness <= std::numeric_limits<scalar_type>::epsilon() ||
            m_material.Z() == 0 || m_material.mass_density() == 0) {
            return false;
        }
        return true;
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
    /// Return the path segment
    DETRAY_HOST_DEVICE
    scalar_type path_segment(const line_plane_intersection& is) const {
        return m_thickness / is.cos_incidence_angle;
    }
    /// Return the path segment in X0
    DETRAY_HOST_DEVICE
    scalar_type path_segment_in_X0(const line_plane_intersection& is) const {
        return m_thickness_in_X0 / is.cos_incidence_angle;
    }
    /// Return the path segment in L0
    DETRAY_HOST_DEVICE
    scalar_type path_segment_in_L0(const line_plane_intersection& is) const {
        return m_thickness_in_L0 / is.cos_incidence_angle;
    }

    private:
    material_type m_material = {};
    scalar_type m_thickness = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_X0 = std::numeric_limits<scalar>::epsilon();
    scalar_type m_thickness_in_L0 = std::numeric_limits<scalar>::epsilon();
};

}  // namespace detray