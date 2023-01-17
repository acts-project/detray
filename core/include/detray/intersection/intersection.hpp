/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"

// System include(s)
#include <cmath>
#include <limits>
#include <ostream>

namespace detray {

namespace intersection {

/// Intersection direction with respect to the normal of the surface
enum class direction {
    e_undefined = -1,  //!< the undefined direction at intersection
    e_opposite = 0,    //!< opposite the surface normal at the intersection
    e_along = 1        //!< along the surface normal at the intersection
};

/// Intersection status: outside, missed, inside, hit ( w/o maks status)
enum class status {
    e_outside = -1,   //!< surface hit but outside
    e_missed = 0,     //!< surface missed
    e_undefined = 1,  //!< surface hit but status not checked
    e_inside = 2      //!< surface hit and inside confirmed
};

}  // namespace intersection

/// @brief This class holds the intersection information
///
/// @tparam point3 is the type of global intsersection vector
/// @tparam point2 is the type of the local intersection vector
template <typename surface_handle_t,
          typename algebra_t = __plugin::transform3<detray::scalar>>
struct line_plane_intersection {

    using scalar_t = typename algebra_t::scalar_type;
    using point3 = typename algebra_t::point3;
    using point2 = typename algebra_t::point2;

    /// Result of the intersection
    intersection::status status{intersection::status::e_undefined};

    /// Direction of the intersection with respect to the track
    intersection::direction direction{intersection::direction::e_undefined};

    /// Distance between track and candidate
    scalar_t path{std::numeric_limits<scalar_t>::infinity()};

    /// cosine of incidence angle
    scalar_t cos_incidence_angle{std::numeric_limits<scalar_t>::infinity()};

    /// Navigation information (next volume to go to)
    dindex volume_link{dindex_invalid};

    /// Handle of the surface this intersection belongs to
    surface_handle_t surface{};

    /// Intersection point in global 3D cartesian coordinate frame
    point3 p3{std::numeric_limits<scalar_t>::infinity(),
              std::numeric_limits<scalar_t>::infinity(),
              std::numeric_limits<scalar_t>::infinity()};

    /// Local position of the intersection on the surface
    point2 p2{std::numeric_limits<scalar_t>::infinity(),
              std::numeric_limits<scalar_t>::infinity()};

    /// @param rhs is the right hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool operator<(const line_plane_intersection &rhs) const {
        return (std::abs(path) < std::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool operator>(const line_plane_intersection &rhs) const {
        return (std::abs(path) > std::abs(rhs.path));
    }

    /// @param rhs is the left hand side intersection for comparison
    DETRAY_HOST_DEVICE
    bool operator==(const line_plane_intersection &rhs) const {
        return std::abs(path - rhs.path) <
               std::numeric_limits<scalar_t>::epsilon();
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    friend std::ostream &operator<<(std::ostream &out_stream,
                                    const line_plane_intersection &is) {
        scalar_t r{getter::perp(is.p3)};
        out_stream << "dist:" << is.path << " [glob: r:" << r
                   << ", z:" << is.p3[2] << " | loc: " << is.p2[0] << ", "
                   << is.p2[1] << "], (sf index:" << is.surface.barcode()
                   << ", links to vol:" << is.volume_link << ")";
        switch (is.status) {
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
        };
        out_stream << std::endl;
        return out_stream;
    }
};

}  // namespace detray
