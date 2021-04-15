/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <climits>
#include <tuple>
#include <optional>

namespace detray
{

    /** Intersection direction with respect to the
     * normal of the surface
     */
    enum intersectiondirection : int
    {
        e_undefined = -1, //!< the undefined direction at intersection
        e_opposite = 0,   //!< opposite the surface normal at the intersection
        e_along = 1       //!< along the surface normal at the intersection
    };

    /** Intersection status: outside, missed, insidem, hit ( w/o maks status) */
    enum intersection_status : int
    {
        e_outside = -1, //!< surface hit but outside
        e_missed = 0,   //!< surface missed
        e_hit = 1,      //!< surface hit but status not checked
        e_inside = 2    //!< surface hit and inside confirmed
    };

    /** This templated class holds the intersection information
     * 
     * @tparam point3_type is the type of global intsersection vector
     * @tparam point2_type is the type of the local intersection vector
     * 
     **/
    template <typename point3_type, typename point2_type>
    struct intersection
    {

        scalar path = std::numeric_limits<scalar>::infinity();
        point3_type point3 = point3_type{std::numeric_limits<scalar>::infinity(),
                                         std::numeric_limits<scalar>::infinity(),
                                         std::numeric_limits<scalar>::infinity()};

        std::optional<point2_type> point2 = std::nullopt;
        intersection_status status = e_missed;
        intersectiondirection direction = e_undefined;
        int index = -1;

        /** @param rhs is the right hand side intersection for comparison 
         **/
        bool operator<(
            const intersection<point3_type, point2_type> &rhs) const
        {
            return (path < rhs.path);
        }

        /** @param rhs is the left hand side intersection for comparison 
         **/
        bool operator>(
            const intersection<point3_type, point2_type> &rhs) const
        {
            return (path > rhs.path);
        }
    };

} // namespace detray