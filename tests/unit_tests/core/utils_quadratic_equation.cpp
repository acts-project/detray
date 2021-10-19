/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

// detray test
#include "tests/common/test_defs.hpp"

// detray core
#include "utils/quadratic_equation.hpp"

using namespace detray;

// This tests the convenience quadratic equation struct
TEST(utils, quad_equation) {
    quadratic_equation<scalar> qe = {{2., 5., -3.}};
    auto solution = qe();

    ASSERT_EQ(std::get<0>(solution), 2);
    ASSERT_NEAR(std::get<1>(solution)[0], -3., 1e-5);
    ASSERT_NEAR(std::get<1>(solution)[1], 0.5, 1e-5);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
