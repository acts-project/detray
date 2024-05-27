/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/detail/algebra.hpp"
#include "detray/definitions/detail/boolean.hpp"
#include "detray/definitions/detail/indexing.hpp"
#include "detray/definitions/detail/math.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <cstdint>
#include <limits>
#include <ostream>

namespace detray {

/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t, typename algebra_t>
struct intersection2D {

    using T = typename algebra_t::value_type;
    using algebra_type = algebra_t;
    using bool_t = dbool<algebra_t>;
    using scalar_type = dscalar<algebra_t>;
    using point2_type = dpoint2D<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using nav_link_t = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t sf_desc{};

    /// Local position of the intersection on the surface
    point3_type local{detail::invalid_value<T>(), detail::invalid_value<T>(),
                      detail::invalid_value<T>()};

    /// Distance between track and candidate
    scalar_type path = detail::invalid_value<T>();

    /// Cosine of incidence angle
    scalar_type cos_incidence_angle = detail::invalid_value<T>();

    /// Navigation information (next volume to go to)
    nav_link_t volume_link{detail::invalid_value<nav_link_t>()};

    /// Result of the intersection (true = inside, false = outside)
    bool_t status{false};

    /// Direction of the intersection with respect to the track (true = along,
    /// false = opposite)
    bool_t direction{true};

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator<(const intersection2D &rhs) const {
        return (math::fabs(path) < math::fabs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator>(const intersection2D &rhs) const {
        return (math::fabs(path) > math::fabs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator==(const intersection2D &rhs) const {
        return math::fabs(path - rhs.path) <
               std::numeric_limits<float>::epsilon();
    }

    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const { return detail::any_of(this->status); }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path
                   << "\tsurface: " << is.sf_desc.barcode()
                   << ", type: " << static_cast<int>(is.sf_desc.mask().id())
                   << ", links to vol:" << is.volume_link << ")"
                   << ", loc [" << is.local[0] << ", " << is.local[1] << ", "
                   << is.local[2] << "], ";

        if constexpr (std::is_scalar_v<bool_t>) {
            out_stream << (is.status ? ", status: inside"
                                     : ", status: outside");
            out_stream << (is.direction ? ", status: along"
                                        : ", status: opposite");
        } else {
            out_stream << ", status: " << is.status;
            out_stream << ", direction: " << is.direction;
        }

        out_stream << std::endl;
        return out_stream;
    }
};

}  // namespace detray
