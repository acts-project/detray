/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cmath>

#include "detray/definitions/qualifiers.hpp"
#include "detray/definitions/invalid_values.hpp"
#include "detray/utils/containers.hpp"

#ifdef DETRAY_CUSTOM_SCALARTYPE
using detray_scalar = DETRAY_CUSTOM_SCALARTYPE;
#else
using detray_scalar = double;
#endif

#define __plugin test

namespace detray {
using scalar = detray_scalar;

namespace test {
using point2 = darray<scalar, 2>;
using point3 = darray<scalar, 3>;

}  // namespace test

namespace getter {

/** Define the perpendicular length
 * @param  is the input vector
 * @return a scalar type */
template <typename point_type>
scalar perp(const point_type &p) {
    return std::sqrt(p[0] * p[0] + p[1] * p[1]);
}

}  // namespace getter

inline test::point2 operator-(const test::point2 &a, const test::point2 &b) {
    return {a[0] - b[0], a[1] - b[1]};
}

DETRAY_HOST_DEVICE
inline bool operator==(const test::point3 &lhs, const test::point3 &rhs) {
    for (int i = 0; i < 3; i++) {
        if (lhs[i] != rhs[i])
            return false;
    }
    return true;
}

// invalid value specialization for test::point3
template <>
DETRAY_HOST_DEVICE inline test::point3 detray::invalid_value() {
    return test::point3{0, 0, 0};
}

}  // namespace detray
