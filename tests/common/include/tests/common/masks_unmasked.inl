/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/masks.hpp"

using namespace detray;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of an unmasked plane
TEST(mask, unmasked) {
    typename mask<unmasked>::loc_point_t p2 = {0.5f, -9.f};

    mask<unmasked> u{};

    ASSERT_TRUE(u.is_inside(p2, 0.f) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = u.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }
}
