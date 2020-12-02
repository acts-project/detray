/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <array>
#include <cmath>

namespace detray
{

    using scalar = float;
    using point2 = std::array<float, 2>;
    using point3 = std::array<float, 3>;

    namespace vector {
        
        /** Define the perpendicular length 
         * @param  is the input vector
         * @return a scalar type */
        scalar perp(const point2& p) {
            return std::sqrt(p[0]*p[0]+p[1]*p[1]);
        }

    }

} // namespace detray
