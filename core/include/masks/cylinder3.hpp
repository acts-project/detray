/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */
#pragma once

#include "core/intersection.hpp"
#include "utils/containers.hpp"
#include "tools/cylinder_intersector.hpp"

#include <cmath>
#include <climits>
#include <optional>

namespace detray
{
    /** This is a simple mask for a full cylinder
     * 
     * @tparam scalar_type the primitive scalar type
     * @tparam kRadialCheck is a boolean to steer wheter the radius compatibility needs to be checked
     * 
     * It is defined by r and the half length.
     **/
    template <typename scalar_type,
              bool kRadialCheck = true,
              typename intersector_type = detray::cylinder_intersector,
              typename links_type = bool,
              unsigned int kMaskIdentifier = 3>
    struct cylinder3
    {
        darray<scalar_type, 3> _values =
            {std::numeric_limits<scalar_type>::infinity(),
             -std::numeric_limits<scalar_type>::infinity(),
             std::numeric_limits<scalar_type>::infinity()};

        links_type _links;

        static constexpr unsigned int mask_identifier = kMaskIdentifier;

        /** Assignment operator from an array, convenience function
         * 
         * @param rhs is the right hand side object
         **/
        cylinder3<scalar_type, kRadialCheck, intersector_type, links_type, kMaskIdentifier> &
        operator=(const darray<scalar_type, 3> &rhs)
        {
            _values = rhs;
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
                                       scalar_type t0 = std::numeric_limits<scalar_type>::epsilon(),
                                       scalar_type t1 = std::numeric_limits<scalar_type>::epsilon()) const
        {
            if (kRadialCheck)
            {
                scalar_type r = getter::perp(p);
                if (std::abs(r - _values[0]) >= t0 + 5 * std::numeric_limits<scalar_type>::epsilon())
                {
                    return e_missed;
                }
            }
            return (_values[1] - t1 <= p[2] and p[2] <= _values[2] + t1) ? e_inside : e_outside;
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
        bool operator==(const cylinder3<scalar_type> &rhs)
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
        intersector_type intersector() { return intersector_type{}; };

        /** Return the volume link */
        const links_type &links() const { return _links; }
    };

} // namespace detray
