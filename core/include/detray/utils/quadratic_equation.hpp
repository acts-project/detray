/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <algorithm>
#include <climits>
#include <cmath>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Struct to solve a quadratic equation of type p[0] * x^2 + p[1] * x + p[2] =
 * 0
 */
template <typename scalar_t,
          template <typename, std::size_t> class array_t = darray,
          template <typename...> class tuple_t = dtuple>
struct quadratic_equation {
    array_t<scalar_t, 3> _params = {0.f, 0.f, 0.f};

    /** Solve the quadratic equation
     **/
    DETRAY_HOST_DEVICE
    tuple_t<int, array_t<scalar_t, 2>> operator()() const {
        scalar_t discriminant =
            _params[1] * _params[1] - 4.f * _params[0] * _params[2];
        if (discriminant < 0.f) {
            return {0,
                    {std::numeric_limits<scalar_t>::infinity(),
                     std::numeric_limits<scalar_t>::infinity()}};
        } else {
            int solutions = (discriminant == 0.f) ? 1 : 2;
            scalar_t q =
                -0.5f *
                (_params[1] + (_params[1] > 0.f ? std::sqrt(discriminant)
                                                : -std::sqrt(discriminant)));
            scalar_t first = q / _params[0];
            scalar_t second = _params[2] / q;
            array_t<scalar_t, 2> poles =
                first < second ? array_t<scalar_t, 2>{first, second}
                               : array_t<scalar_t, 2>{second, first};
            return {solutions, poles};
        }
    }
};

}  // namespace detray
