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

namespace detray {

namespace intersection {

/// Intersection direction with respect to the normal of the surface
enum class direction : std::uint_least8_t {
    e_undefined = 0u,  //!< the undefined direction at intersection
    e_opposite = 1u,   //!< opposite the surface normal at the intersection
    e_along = 2u       //!< along the surface normal at the intersection
};

/// Intersection status: outside, missed, inside, hit ( w/o maks status)
enum class status : std::uint_least8_t {
    e_outside = 0u,    //!< surface hit but outside
    e_missed = 1u,     //!< surface missed
    e_undefined = 2u,  //!< surface hit but status not checked
    e_inside = 3u      //!< surface hit and inside confirmed
};

}  // namespace intersection

#if (IS_SOA)
/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t, typename T = scalar,
          template <typename> class algebra_t = vc_soa>
struct intersection2D {

    using scalar_t = dscalar<algebra_t<T>>;
    using bool_mask = dbool<algebra_t<T>>;
    using point3D = dpoint3D<algebra_t<T>>;
    using vector3D = dvector3D<algebra_t<T>>;
    using transform3D = dtransform3D<algebra_t<T>>;
    using nav_link_t = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t surface;

    /// Local position of the intersection on the surface
    point3D local{detail::invalid_value<T>(), detail::invalid_value<T>(),
                  detail::invalid_value<T>()};

    /// Distance between track and candidate
    scalar_t path = detail::invalid_value<T>();

    /// Cosine of incidence angle
    scalar_t cos_incidence_angle = detail::invalid_value<T>();

    /// Navigation information (next volume to go to)
    nav_link_t volume_link = detail::invalid_value<nav_link_t>();

    /// Result of the intersection (true = inside, false = outside)
    bool_mask status = false;

    /// Direction of the intersection with respect to the track (true = along,
    /// false = opposite)
    bool_mask direction = false;

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    auto operator<(const intersection2D &rhs) const {
        return (Vc::abs(path) < Vc::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    auto operator>(const intersection2D &rhs) const {
        return (Vc::abs(path) > Vc::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    auto operator==(const intersection2D &rhs) const {
        return Vc::abs(path - rhs.path) < std::numeric_limits<T>::epsilon();
    }

    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const { return Vc::any_of(this->status); }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path << ", loc pos [" << is.local[0] << ", "
                   << is.local[1] << ", " << is.local[2] << "]"
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
#else
/// @brief This class holds the intersection information.
///
/// @tparam surface_descr_t is the type of surface descriptor
template <typename surface_descr_t,
          typename algebra_t = __plugin::transform3<detray::scalar>>
struct intersection2D {

    using transform3_type = algebra_t;
    using scalar_type = typename algebra_t::scalar_type;
    using point3 = typename algebra_t::point3;
    using point2 = typename algebra_t::point2;
    using nav_link_type = typename surface_descr_t::navigation_link;

    /// Descriptor of the surface this intersection belongs to
    surface_descr_t sf_desc;

    /// Local position of the intersection on the surface
    point3 local{detail::invalid_value<scalar_type>(),
                 detail::invalid_value<scalar_type>(),
                 detail::invalid_value<scalar_type>()};

    /// Distance between track and candidate
    scalar_type path{detail::invalid_value<scalar_type>()};

    /// cosine of incidence angle
    scalar_type cos_incidence_angle{detail::invalid_value<scalar_type>()};

    /// Navigation information (next volume to go to)
    nav_link_type volume_link{detail::invalid_value<nav_link_type>()};

    /// Result of the intersection
    intersection::status status{intersection::status::e_undefined};

    /// Direction of the intersection with respect to the track
    intersection::direction direction{intersection::direction::e_undefined};

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
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
    }

    DETRAY_HOST_DEVICE
    constexpr bool is_inside() const {
        return this->status == intersection::status::e_inside;
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const intersection2D &is) {
        out_stream << "dist:" << is.path
                   << "\tsurface: " << is.sf_desc.barcode()
                   << ", links to vol:" << is.volume_link << ")"
                   << ", loc [" << is.local[0] << ", " << is.local[1] << ", "
                   << is.local[2] << "], ";

        switch (is.status) {
            case intersection::status::e_outside:
                out_stream << "status: outside";
                break;
            case intersection::status::e_missed:
                out_stream << "status: missed";
                break;
            case intersection::status::e_undefined:
                out_stream << "status: undefined";
                break;
            case intersection::status::e_inside:
                out_stream << "status: inside";
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
        };

        return out_stream;
    }
};
#endif

}  // namespace detray
