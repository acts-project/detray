/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Licenced under: Apache-2, see LICENSE file
 */
#pragma once

#include <climits>
#include <optional>
#include <tuple>

namespace detray
{

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
     * @tparam scalar_type is the type of the used scalar for intersecting
     * @tparam point3_type is the type of global intsersection vector
     * @tparam point2_type is the type of the local intersection vector
     * 
     **/
    template <typename scalar_type, typename point3_type, typename point2_type>
    struct intersection
    {

        scalar_type _path = std::numeric_limits<scalar_type>::infinity();
        point3_type _point3 = point3_type(std::numeric_limits<scalar_type>::infinity(), 
                                          std::numeric_limits<scalar_type>::infinity(), 
                                          std::numeric_limits<scalar_type>::infinity());
        std::optional<point2_type> _point2 = std::nullopt;
        intersection_status _status = e_missed;

        /** @param rhs is the right hand side intersection for comparison 
         **/
        bool operator<(
            const intersection<scalar_type, point3_type, point2_type> &rhs) const
        {
            return (_path < rhs._path);
        }

        /** @param rhs is the left hand side intersection for comparison 
         **/
        bool operator>(
            const intersection<scalar_type, point3_type, point2_type> &rhs) const
        {
            return (_path > rhs._path);
        }
    };

} // namespace detray