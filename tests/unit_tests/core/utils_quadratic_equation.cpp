/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "detray/definitions/algebra.hpp"
#include "detray/utils/quadratic_equation.hpp"

using namespace detray;

// This tests the convenience quadratic equation struct
TEST(utils, quad_equation) {
    quadratic_equation<scalar> qe = {{2.f, 5.f, -3.f}};
    auto solution = qe();

    ASSERT_EQ(std::get<0>(solution), 2);
    ASSERT_NEAR(std::get<1>(solution)[0], -3.f, 1e-7f);
    ASSERT_NEAR(std::get<1>(solution)[1], 0.5f, 1e-7f);
}
