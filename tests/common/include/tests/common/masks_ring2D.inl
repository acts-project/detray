/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;

/// This tests the basic functionality of a ring
TEST(mask, ring2D) {
    using point_t = typename mask<ring2D<>>::loc_point_t;

    point_t p2_pl_in = {0.5, -2.};
    point_t p2_pl_edge = {0., 3.5};
    point_t p2_pl_out = {3.6, 5.};

    constexpr scalar inner_r{0. * unit_constants::mm};
    constexpr scalar outer_r{3.5 * unit_constants::mm};

    mask<ring2D<>> r2{0UL, inner_r, outer_r};

    ASSERT_FLOAT_EQ(r2[ring2D<>::e_inner_r], static_cast<scalar>(0.));
    ASSERT_FLOAT_EQ(r2[ring2D<>::e_outer_r], static_cast<scalar>(3.5));

    ASSERT_TRUE(r2.is_inside(p2_pl_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_pl_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_pl_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_pl_out, 1.2) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = r2.projection_matrix();
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(r2)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}