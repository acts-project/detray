/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/geometry/barcode.hpp"

// System include(s)
#include <climits>
#include <cmath>
#include <sstream>

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
    e_outside = -1,  //!< surface hit but outside
    e_missed = 0,    //!< surface missed
    e_hit = 1,       //!< surface hit but status not checked
    e_inside = 2     //!< surface hit and inside confirmed
};

}  // namespace intersection

/// @brief This class holds the intersection information
///
/// @tparam point3 is the type of global intsersection vector
/// @tparam point2 is the type of the local intersection vector
struct line_plane_intersection {

    using point3 = __plugin::point3<scalar>;
    using vector3 = __plugin::vector3<scalar>;
    using point2 = __plugin::point2<scalar>;

    /// Distance between track and candidate
    scalar path{std::numeric_limits<scalar>::infinity()};
    point3 p3{std::numeric_limits<scalar>::infinity(),
              std::numeric_limits<scalar>::infinity(),
              std::numeric_limits<scalar>::infinity()};

    /// Local position of the intersection on the surface
    point2 p2{std::numeric_limits<scalar>::infinity(),
              std::numeric_limits<scalar>::infinity()};

    /// Result of the intersection
    intersection::status status{intersection::status::e_missed};

    /// Direction of the intersection with respect to the track
    intersection::direction direction{intersection::direction::e_undefined};

    /// Mask index
    dindex mask_index{dindex_invalid};

    /// Primitive this intersection belongs to
    geometry::barcode barcode{};

    /// Navigation information
    dindex volume_link{dindex_invalid};

    /// Surface id for this intersection (sensitive, portal, passive)
    surface_id sf_id{surface_id::e_sensitive};

    // cosine of incidence angle
    scalar cos_incidence_angle{1.f};

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
               std::numeric_limits<scalar>::epsilon();
    }

    /// Transform to a string for output debugging
    DETRAY_HOST
    std::string to_string() const {
        std::stringstream out_stream;
        scalar r{getter::perp(p3)};
        out_stream << "dist:" << path << " [r:" << r << ", z:" << p3[2]
                   << "], (sf index:" << barcode
                   << ", links to vol:" << volume_link << ")";
        switch (status) {
            case intersection::status::e_outside:
                out_stream << ", status: outside";
                break;
            case intersection::status::e_missed:
                out_stream << ", status: missed";
                break;
            case intersection::status::e_hit:
                out_stream << ", status: hit";
                break;
            case intersection::status::e_inside:
                out_stream << ", status: inside";
                break;
        };
        switch (direction) {
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
        return out_stream.str();
    }
};

}  // namespace detray
