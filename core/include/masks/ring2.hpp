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
    /** This is a simple 2-dimensional mask for a closed ring
     * 
     * It is defined by the two radii _v[0] and  _v[1], 
     * and can be checked with a tolerance in t0 and t1.
     **/
    template <typename scalar_type, 
              typename intersector_type = planar_intersector, 
              typename links_type = bool,
              unsigned int kMaskIdentifier=2>
    struct ring2
    {
        darray<scalar_type, 2> _v =
            {0.,
             std::numeric_limits<scalar_type>::infinity()};

        links_type _links;

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        ring2<scalar_type, intersector_type, links_type, kMaskIdentifier>&
        operator=(const darray<scalar_type, 2> &rhs)
        {
            _v = rhs;
            return (*this);
        }

        /** Mask operation 
         * 
         * @tparam point2_type is the type of the point to be checked w.r.t. to
         * the mask bounds
         * 
         * @param p the point to be checked
         * @param t0 is the tolerance in local 0
         * @param t1 is the tolerance in local 1 and is ignored
         * 
         * @return an intersection status e_inside / e_outside
         **/
        template <typename point2_type>
        intersection_status operator()(const point2_type &p,
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
            scalar_type r = getter::perp(p - point2_type{t0, t1});
            return (r + t0 >= _v[0] and r <= _v[1] + t0) ? e_inside : e_outside;
        }

        /** Equality operator from an array, convenience function
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const darray<scalar_type, 2> &rhs)
        {
            return (std::abs(_v[0] - rhs[0]) < std::numeric_limits<scalar_type>::epsilon() and std::abs(_v[1] - rhs[1]) < std::numeric_limits<scalar_type>::epsilon());
        }

        /** Equality operator 
         * 
         * @param rhs is the rectangle to be compared with
         * 
         * checks identity within epsilon and @return s a boolean*
         **/
        bool operator==(const ring2<scalar_type> &rhs)
        {
            return operator==(rhs._v);
        }

        /** Access operator - non-const
         * @return the reference to the member variable
         */
        scalar_type& operator[](unsigned int value_index) {
            return _v[value_index];
        }

        /** Access operator - non-const
         * @return a copy of the member variable
         */
        scalar_type operator[](unsigned int value_index) const {
            return _v[value_index];
        }

        /** Return an associated intersector type */
        intersector_type intersector() { return intersector_type{}; };

        /** Mask identifier */
        constexpr unsigned int mask_identifier() { return kMaskIdentifier; }

        /** Return the volume link */
        const links_type& links() const { return _links; }
    };

} // namespace detray
