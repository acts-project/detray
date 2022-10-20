/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/axis_rotation.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;
using vector3 = __plugin::vector3<scalar>;
using transform3 = __plugin::transform3<detray::scalar>;
using matrix_operator = typename transform3::matrix_actor;
using size_type = typename matrix_operator::size_ty;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

const scalar isclose = 1e-5;

TEST(utils, axis_rotation) {

    const vector3 axis{0, 0, 1};

    const vector3 v1{1, 0, 0};

    const auto u1 = axis_rotation<transform3>(axis, M_PI_2)(v1);

    EXPECT_NEAR(u1[0], 0, isclose);
    EXPECT_NEAR(u1[1], 1, isclose);
    EXPECT_NEAR(u1[2], 0, isclose);

    matrix_type<3, 1> v2;
    matrix_operator().element(v2, 0, 0) = 1. / M_SQRT2;
    matrix_operator().element(v2, 1, 0) = 1. / M_SQRT2;
    matrix_operator().element(v2, 2, 0) = 0;

    const auto u2 = axis_rotation<transform3>(axis, M_PI_4)(v2);

    EXPECT_NEAR(matrix_operator().element(u2, 0, 0), 0, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 1, 0), 1, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 2, 0), 0, isclose);
}
