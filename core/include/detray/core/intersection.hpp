/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <optional>
#include <tuple>

#include "detray/utils/indexing.hpp"

namespace detray {

using point3 = __plugin::point3;
using vector3 = __plugin::vector3;
using point2 = __plugin::point2;

/** Intersection direction with respect to the
 * normal of the surface
 */
enum intersection_direction : int {
    e_undefined = -1,  //!< the undefined direction at intersection
    e_opposite = 0,    //!< opposite the surface normal at the intersection
    e_along = 1        //!< along the surface normal at the intersection
};

/** Intersection status: outside, missed, insidem, hit ( w/o maks status) */
enum intersection_status : int {
    e_outside = -1,  //!< surface hit but outside
    e_missed = 0,    //!< surface missed
    e_hit = 1,       //!< surface hit but status not checked
    e_inside = 2     //!< surface hit and inside confirmed
};

/** This templated class holds the intersection information
 *
 * @tparam point3 is the type of global intsersection vector
 * @tparam point2 is the type of the local intersection vector
 *
 **/
struct intersection {

    scalar path = std::numeric_limits<scalar>::infinity();
    point3 p3 = point3{std::numeric_limits<scalar>::infinity(),
                       std::numeric_limits<scalar>::infinity(),
                       std::numeric_limits<scalar>::infinity()};

    point2 p2 = {std::numeric_limits<scalar>::infinity(),
                 std::numeric_limits<scalar>::infinity()};
    intersection_status status = e_missed;
    intersection_direction direction = e_undefined;
    dindex index = dindex_invalid;
    dindex link = dindex_invalid;

    /** @param rhs is the right hand side intersection for comparison
     **/
    bool operator<(const intersection &rhs) const { return (path < rhs.path); }

    /** @param rhs is the left hand side intersection for comparison
     **/
    bool operator>(const intersection &rhs) const { return (path > rhs.path); }

    /** @param rhs is the left hand side intersection for comparison
     **/
    bool operator==(const intersection &rhs) const {
        return (path == rhs.path);
    }
};

}  // namespace detray
