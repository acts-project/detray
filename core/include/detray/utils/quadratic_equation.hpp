/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <tuple>

#include "detray/definitions/qualifiers.hpp"

namespace detray {

/** Struct to solve a quadratic equation of type p[0] * x^2 + p[1] * x + p[2] =
 * 0
 */
template <typename scalar_type,
          template <typename, unsigned int> class array_type = darray,
          template <typename...> class tuple_type = dtuple>
struct quadratic_equation {
    array_type<scalar_type, 3> _params = {0., 0., 0.};

    /** Solve the quadratic equation
     **/
    DETRAY_HOST_DEVICE
    tuple_type<int, array_type<scalar_type, 2>> operator()() const {
        scalar_type discriminant =
            _params[1] * _params[1] - 4 * _params[0] * _params[2];
        if (discriminant < 0.) {
            return {0,
                    {std::numeric_limits<scalar_type>::infinity(),
                     std::numeric_limits<scalar_type>::infinity()}};
        } else {
            int solutions = (discriminant == 0.) ? 1 : 2;
            double q = -0.5 * (_params[1] + (_params[1] > 0
                                                 ? std::sqrt(discriminant)
                                                 : -std::sqrt(discriminant)));
            scalar_type first = q / _params[0];
            scalar_type second = _params[2] / q;
            array_type<scalar_type, 2> poles =
                first < second ? array_type<scalar_type, 2>{first, second}
                               : array_type<scalar_type, 2>{second, first};
            return {solutions, poles};
        }
    }
};

}  // namespace detray
