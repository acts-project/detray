/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray {

/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t, typename T,
          template <typename> class algebra_t>
struct intersection2D {

    using value_t = T;
    using scalar_t = dscalar<algebra_t<T>>;
    using bool_t = dbool<algebra_t<T>>;
    using point2D = dpoint2D<algebra_t<T>>;
    using point3D = dpoint3D<algebra_t<T>>;
    using vector3D = dvector3D<algebra_t<T>>;
    using transform3D = dtransform3D<algebra_t<T>>;
    using nav_link_t = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t sf_desc;

    /// Local position of the intersection on the surface
    point3D local{detail::invalid_value<T>(), detail::invalid_value<T>(),
                  detail::invalid_value<T>()};

    /// Distance between track and candidate
    scalar_t path = detail::invalid_value<T>();

    /// Cosine of incidence angle
    scalar_t cos_incidence_angle = detail::invalid_value<T>();

    /// Navigation information (next volume to go to)
    nav_link_t volume_link;

    /// Result of the intersection (true = inside, false = outside)
    bool_t status;

    /// Direction of the intersection with respect to the track (true = along,
    /// false = opposite)
    bool_t direction;

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator<(const intersection2D &rhs) const {
        return (math_ns::abs(path) < math_ns::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator>(const intersection2D &rhs) const {
        return (math_ns::abs(path) > math_ns::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool_t operator==(const intersection2D &rhs) const {
        return math_ns::abs(path - rhs.path) <
               std::numeric_limits<T>::epsilon();
    }

    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const { return detail::any_of(this->status); }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path << ", loc pos [" << is.local[0] << ", "
                   << is.local[1] << ", " << is.local[2] << "]"
                   << ", (sf index:" << is.sf_desc.barcode()
                   << ", links to vol:" << is.volume_link << ")";

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
