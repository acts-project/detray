/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/utils/axis_rotation.hpp"

#include "detray/definitions/units.hpp"
#include "detray/test/types.hpp"

// GTest include(s).
#include <gtest/gtest.h>

using namespace detray;

using matrix_operator = typename test::transform3::matrix_actor;
using size_type = typename matrix_operator::size_ty;
template <size_type ROWS, size_type COLS>
using matrix_type = typename matrix_operator::template matrix_type<ROWS, COLS>;

static constexpr scalar isclose{1e-5f};

GTEST_TEST(detray_utils, axis_rotation) {

    const test::vector3 axis{0.f, 0.f, 3.f};

    const test::vector3 v1{1.f, 0.f, 0.f};

    const auto u1 =
        axis_rotation<test::transform3>(axis, constant<scalar>::pi_2)(v1);

    EXPECT_NEAR(u1[0], 0.f, isclose);
    EXPECT_NEAR(u1[1], 1.f, isclose);
    EXPECT_NEAR(u1[2], 0.f, isclose);

    matrix_type<3, 1> v2;
    matrix_operator().element(v2, 0, 0) = constant<scalar>::inv_sqrt2;
    matrix_operator().element(v2, 1, 0) = constant<scalar>::inv_sqrt2;
    matrix_operator().element(v2, 2, 0) = 0.f;

    const auto u2 =
        axis_rotation<test::transform3>(axis, constant<scalar>::pi_4)(v2);

    EXPECT_NEAR(matrix_operator().element(u2, 0, 0), 0.f, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 1, 0), 1.f, isclose);
    EXPECT_NEAR(matrix_operator().element(u2, 2, 0), 0.f, isclose);
}
