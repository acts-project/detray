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
    /** This is a simple mask for single parameter bound mask
     * 
     **/
    template <typename scalar_type, 
              unsigned int kCheckIndex, 
              typename intersector_type = planar_intersector, 
              typename links_type = bool,
              unsigned int kMaskIdentifier=4>
    struct single3
    {

        using mask_values = darray<scalar_type, 1>;

        mask_values _values=
            {std::numeric_limits<scalar_type>::infinity()};

        links_type _links;

        static constexpr unsigned int mask_identifier = kMaskIdentifier;

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        single3<scalar_type, kCheckIndex, intersector_type, links_type, kMaskIdentifier>&
        operator=(const darray<scalar_type, 1> &rhs)
        {
            _values= rhs;
            return (*this);
        }

        /** Mask operation 
         * 
         * @tparam point3_type is the type of the point to be checked w.r.t. to
         * the mask bounds, it's assumed to be within the cylinder 3D frame
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0 (radius)
         * @param t1 is the tolerance in local 1 (z)
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point3_type>
        intersection_status operator()(const point3_type &p,
                                       scalar_type t = std::numeric_limits<scalar_type>::epsilon()) const
        {     
            return (std::abs(p[kCheckIndex]) <= _values[0] + t) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (_values== rhs);
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const single3<scalar_type, kCheckIndex> &rhs)
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
        const links_type& links() const { return _links; }
    };

} // namespace detray
