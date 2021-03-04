/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "tools/planar_intersector.hpp"

#include <cmath>
#include <climits>

namespace detray
{

    /** This is a simple 2-dimensional mask for a regular rectangle
     * 
     * It is defined by half length in local0 coordinates _values[0] and _values[1], 
     * and can be checked with a tolerance in t0 and t1.
     **/
    template <typename scalar_type,
              typename intersector_type = planar_intersector,
              typename links_type = bool,
              unsigned int kMaskIdentifier = 0>
    struct rectangle2
    {
        using mask_tolerance = darray<scalar_type, 2>;

        using mask_values = darray<scalar_type, 2>;

        mask_values _values =
            {std::numeric_limits<scalar_type>::infinity(),
             std::numeric_limits<scalar_type>::infinity()};

        links_type _links;

        static constexpr unsigned int mask_identifier = kMaskIdentifier;

        static constexpr mask_tolerance within_epsilon = {std::numeric_limits<scalar_type>::epsilon(),
                                                           std::numeric_limits<scalar_type>::epsilon()};

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        rectangle2<scalar_type, intersector_type, links_type, kMaskIdentifier> &
        operator=(const darray<scalar_type, 2> &rhs)
        {
            _values = rhs;
            return (*this);
        }

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked
         * @param t is the tolerance tuple in (l0, l1)
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename local_type>
        intersection_status is_inside(const typename local_type::point2 &p,
                                      const mask_tolerance &t = within_epsilon) const
        {
            return (std::abs(p[0]) <= _values[0] + t[0] and std::abs(p[1]) <= _values[1] + t[1]) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (_values == rhs);
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const rectangle2<scalar_type> &rhs)
        {
            return operator==(rhs._values);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar_type &operator[](unsigned int value_index)
        {
            return _values[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar_type operator[](unsigned int value_index) const
        {
            return _values[value_index];
        }

        /** Return an associated intersector type */
        intersector_type intersector() const { return intersector_type{}; };

        /** Return the volume link - const reference */
        const links_type &links() const { return _links; }

        /** Return the volume link - non-const access */
        links_type &links() { return _links; }
    };

} // namespace detray
