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
#include <iostream>

namespace detray
{
    /** This is a simple 2-dimensional mask for a regular trapezoid
     * 
     * It is defined by half lengths in local0 coordinate _values[0] and _values[1] at
     * -/+ half length in the local1 coordinate _values[2]
     **/
    template <typename scalar_type,
              typename intersector_type = planar_intersector,
              typename links_type = bool,
              unsigned int kMaskIdentifier = 1>
    struct trapezoid2
    {
        using mask_values = darray<scalar_type, 3>;

        mask_values _values =
            {std::numeric_limits<scalar_type>::infinity(),
             std::numeric_limits<scalar_type>::infinity(),
             std::numeric_limits<scalar_type>::infinity()};

        links_type _links;

        static constexpr unsigned int mask_identifier = kMaskIdentifier;

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        trapezoid2<scalar_type, intersector_type, links_type, kMaskIdentifier> &
        operator=(const darray<scalar_type, 3> &rhs)
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
         * @param t0 is the tolerance in local 0
         * @param t1 is the tolerance in local 1
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point2_type>
        intersection_status operator()(const point2_type &p,
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
            scalar_type rel_y = (_values[2] + p[1]) / (2 * _values[2]);
            return (std::abs(p[0]) <= _values[0] + rel_y * (_values[1] - _values[0]) + t0 and std::abs(p[1]) <= _values[2] + t1) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 3> &rhs)
        {
            return (_values == rhs);
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const trapezoid2<scalar_type> &rhs)
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

        /** Return the volume link */
        const links_type &links() const { return _links; }
    };

} // namespace detray
