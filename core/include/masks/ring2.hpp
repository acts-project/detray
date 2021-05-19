/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "masks/mask_identifier.hpp"
#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "tools/planar_intersector.hpp"

#include <cmath>
#include <climits>

namespace detray
{
    /** This is a simple 2-dimensional mask for a closed ring
     * 
     * @tparam intersector_type is a struct used for intersecting this cylinder
     * @tparam local_type is the default local frame type
     * @tparam links_type is an object where the mask can link to 
     * @tparam kMaskContext is a unique mask identifier in a certain context
     *  
     * It is defined by the two radii _values[0] and  _values[1], 
     * and can be checked with a tolerance in t[0] and t[1].
     * 
     * @note  While the mask_context can change depending on the typed container
     * structure the mask_identifier is a const expression that determines the
     * mask type once for all.
     * 
     **/
    template <typename intersector_type = planar_intersector,
              typename local_type = __plugin::cartesian2,
              typename links_type = bool,
              unsigned int kMaskContext = e_ring2>
    struct ring2
    {
        using mask_tolerance = scalar;

        using mask_values = darray<scalar, 2>;

        using mask_links_type = links_type;

        mask_values _values = {0., std::numeric_limits<scalar>::infinity()};

        links_type _links;

        static constexpr unsigned int mask_context = kMaskContext;

        static constexpr unsigned int mask_identifier = e_ring2;

        static constexpr mask_tolerance within_epsilon = std::numeric_limits<scalar>::epsilon();

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        ring2<intersector_type, local_type, links_type, kMaskContext> &
        operator=(const darray<scalar, 2> &rhs)
        {
            _values = rhs;
            return (*this);
        }

        /** Mask operation 
         * 
         * @tparam inside_local_type is the local type for inside checking
         * 
         * @param p the point to be checked
         * @param t is the tolerance in r
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename inside_local_type>
        intersection_status is_inside(const point2 &p,
                                      const mask_tolerance &t = within_epsilon) const
        {
            if constexpr (std::is_same_v<inside_local_type, __plugin::cartesian2>)
            {
                scalar r = getter::perp(p);
                return (r + t >= _values[0] and r <= _values[1] + t) ? e_inside : e_outside;
            }

            return (p[0] + t >= _values[0] and p[0] <= _values[1] + t) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar, 2> &rhs)
        {
            return (_values == rhs);
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const ring2<> &rhs)
        {
            return operator==(rhs._values);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar &operator[](unsigned int value_index)
        {
            return _values[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar operator[](unsigned int value_index) const
        {
            return _values[value_index];
        }

        /** Return an associated intersector type */
        intersector_type intersector() const { return intersector_type{}; };

        /** Return the values */
        const mask_values &values() const { return _values; }

        /** Return the local frame type */
        local_type local() const { return local_type{}; }

        /** Return the volume link - const reference */
        const links_type &links() const { return _links; }

        /** Return the volume link - non-const access */
        links_type &links() { return _links; }
    };

} // namespace detray
