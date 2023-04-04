/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cmath>

#include "detray/definitions/containers.hpp"
#include "detray/definitions/qualifiers.hpp"
#include "detray/utils/invalid_values.hpp"

#ifdef DETRAY_CUSTOM_SCALARTYPE
using detray_scalar = DETRAY_CUSTOM_SCALARTYPE;
#else
using detray_scalar = double;
#endif

namespace detray {

namespace test {

template <typename T>
using point2 = darray<T, 2>;
template <typename T>
using point3 = darray<T, 3>;

}  // namespace test

template <typename T>
inline test::point2<T> operator-(const test::point2<T> &a,
                                 const test::point2<T> &b) {
    return {a[0] - b[0], a[1] - b[1]};
}

template <typename T>
DETRAY_HOST_DEVICE inline bool operator==(const test::point3<T> &lhs,
                                          const test::point3<T> &rhs) {
    for (int i = 0; i < 3; i++) {
        if (lhs[i] != rhs[i])
            return false;
    }
    return true;
}

}  // namespace detray
