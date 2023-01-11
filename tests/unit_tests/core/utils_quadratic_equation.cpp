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
#include "detray/utils/quadratic_equation.hpp"

using namespace detray;

// This tests the convenience quadratic equation class
TEST(utils, quadratic_equation) {

    static constexpr scalar epsilon{1e-5};

    // No solution
    detail::quadratic_equation<scalar> qe1{1.5, 0., 1.};

    ASSERT_EQ(qe1.solutions(), 0);

    // One solution
    detail::quadratic_equation<scalar> qe2{1., 0., 0.};

    ASSERT_EQ(qe2.solutions(), 1);
    EXPECT_NEAR(qe2.smaller(), 0.f, epsilon);

    detail::quadratic_equation<scalar> qe3{0., 1., 2.};

    ASSERT_EQ(qe3.solutions(), 1);
    EXPECT_NEAR(qe3.smaller(), -2.f, epsilon);

    // Two solutions
    detail::quadratic_equation<scalar> qe4{2., 5., 3.};

    ASSERT_EQ(qe4.solutions(), 2);
    EXPECT_NEAR(qe4.smaller(), -1.5f, epsilon);
    EXPECT_NEAR(qe4.larger(), -1.f, epsilon);

    detail::quadratic_equation<scalar> qe5{2., 5., -3.};

    ASSERT_EQ(qe5.solutions(), 2);
    EXPECT_NEAR(qe5.smaller(), -3.f, epsilon);
    EXPECT_NEAR(qe5.larger(), 0.5f, epsilon);

    detail::quadratic_equation<scalar> qe6{2., -5., 3.};

    ASSERT_EQ(qe6.solutions(), 2);
    EXPECT_NEAR(qe6.smaller(), 1.f, epsilon);
    EXPECT_NEAR(qe6.larger(), 1.5f, epsilon);

    detail::quadratic_equation<scalar> qe7{2., -5., -3.};

    ASSERT_EQ(qe7.solutions(), 2);
    EXPECT_NEAR(qe7.smaller(), -0.5f, epsilon);
    EXPECT_NEAR(qe7.larger(), 3.f, epsilon);

    detail::quadratic_equation<scalar> qe8{2., -5., 0.};

    ASSERT_EQ(qe8.solutions(), 2);
    EXPECT_NEAR(qe8.smaller(), 0.f, epsilon);
    EXPECT_NEAR(qe8.larger(), 5.f / 2.f, epsilon);
}
