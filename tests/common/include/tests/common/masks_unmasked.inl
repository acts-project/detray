/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/masks.hpp"

using namespace detray;

/// This tests the basic functionality of an unmasked plane
TEST(mask, unmasked) {
    typename mask<unmasked>::loc_point_t p2 = {0.5, -9.};

    mask<unmasked> u{};

    ASSERT_TRUE(u.is_inside(p2, 0) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = u.projection_matrix();
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(u)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}
