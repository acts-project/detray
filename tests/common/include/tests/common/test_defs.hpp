/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <cmath>

#ifdef DETRAY_CUSTOM_SCALARTYPE
using detray_scalar = DETRAY_CUSTOM_SCALARTYPE;
#else
using detray_scalar = double;
#endif

#define plugin test

namespace detray
{
    using scalar = detray_scalar;

    namespace test
    {
        using point2 = std::array<scalar, 2>;
        using point3 = std::array<scalar, 3>;
    } // namespace test

    namespace vector
    {

        /** Define the perpendicular length 
         * @param  is the input vector
         * @return a scalar type */
        template <typename point_type>
        scalar perp(const point_type &p)
        {
            return std::sqrt(p[0] * p[0] + p[1] * p[1]);
        }

    } // namespace vector

    test::point2 operator-(const test::point2 &a, const test::point2 &b)
    {
        return {a[0] - b[0], a[1] - b[1]};
    }

} // namespace detray
