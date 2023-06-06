/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <cstdint>
#include <limits>
#include <ostream>

namespace detray::soa {

#if(IS_SOA)
/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t,
          typename algebra_t = __plugin::transform3<detray::scalar>>
struct intersection2D {

    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using value_type = typename algebra_t::value_type;
    using point3 = typename algebra_t::point3;
    using point2 = typename algebra_t::point2;
    using nav_link_type = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t surface;

    /// Local position of the intersection on the surface
    point3 local{detail::invalid_value<scalar_type>(),
                 detail::invalid_value<scalar_type>(),
                 detail::invalid_value<scalar_type>()};

    /// Distance between track and candidate
    value_type path{detail::invalid_value<scalar_type>()};

    /// cosine of incidence angle
    value_type cos_incidence_angle{detail::invalid_value<scalar_type>()};

    /// Navigation information (next volume to go to)
    nav_link_type volume_link{detail::invalid_value<typename nav_link_type::value_type>()};

    /// Result of the intersection (true = inside, false = outside)
    typename value_type::mask_type status{};

    /// Direction of the intersection with respect to the track (true = along, 
    /// false = opposite)
    typename value_type::mask_type direction{};

    /// @param rhs is the right hand side intersection for comparison
    /*DETRAY_HOST_DEVICE
    bool operator<(const intersection2D &rhs) const {
        return (std::abs(path) < std::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool operator>(const intersection2D &rhs) const {
        return (std::abs(path) > std::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool operator==(const intersection2D &rhs) const {
        return std::abs(path - rhs.path) <
               std::numeric_limits<scalar_type>::epsilon();
    }*/

    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const {
        return Vc::any_of(this->status);
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path
                   << ", loc pos [" << is.local[0] << ", " << is.local[1] 
                   << ", " << is.local[2] << "]"
                   << ", (sf index:" << is.surface.barcode()
                   << ", links to vol:" << is.volume_link << ")";
        /*switch (is.status) {
            case intersection::status::e_outside:
                out_stream << ", status: outside";
                break;
            case intersection::status::e_missed:
                out_stream << ", status: missed";
                break;
            case intersection::status::e_undefined:
                out_stream << ", status: undefined";
                break;
            case intersection::status::e_inside:
                out_stream << ", status: inside";
                break;
        };
        switch (is.direction) {
            case intersection::direction::e_undefined:
                out_stream << ", direction: undefined";
                break;
            case intersection::direction::e_opposite:
                out_stream << ", direction: opposite";
                break;
            case intersection::direction::e_along:
                out_stream << ", direction: along";
                break;
        };*/
        out_stream << std::endl;
        return out_stream;
    }
};
#endif

}  // namespace detray
