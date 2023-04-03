/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
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

constexpr scalar isclose{1e-5f};

TEST(utils, axis_rotation) {

    const vector3 axis{0.f, 0.f, 3.f};

    const vector3 v1{1.f, 0.f, 0.f};

    const auto u1 = axis_rotation<transform3>(axis, constant<scalar>::pi_2)(v1);

    EXPECT_NEAR(u1[0], 0.f, isclose);
    EXPECT_NEAR(u1[1], 1.f, isclose);
    EXPECT_NEAR(u1[2], 0.f, isclose);

    matrix_type<3, 1> v2;
    matrix_operator().element(v2, 0, 0) = constant<scalar>::inv_sqrt2;
    matrix_operator().element(v2, 1, 0) = constant<scalar>::inv_sqrt2;
    matrix_operator().element(v2, 2, 0) = 0.f;

    const auto u2 = axis_rotation<transform3>(axis, constant<scalar>::pi_4)(v2);

    EXPECT_NEAR(matrix_operator().element(u2, 0, 0), 0.f, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 1, 0), 1.f, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 2, 0), 0.f, isclose);
}
