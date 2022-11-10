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
using namespace __plugin;

constexpr scalar r{3. * unit<scalar>::mm};
constexpr scalar hz{4. * unit<scalar>::mm};

/// This tests the basic functionality of a 3D cylinder
TEST(mask, cylinder3D) {
    using point_t = typename mask<cylinder3D>::loc_point_t;

    point_t p3_in = {r, 0., -1.};
    point_t p3_edge = {0., r, hz};
    point_t p3_out = {static_cast<scalar>(r / std::sqrt(2.)),
                      static_cast<scalar>(r / std::sqrt(2.)), 4.5};
    point_t p3_off = {1., 1., -9.};

    // Test radius to be on surface, too
    mask<cylinder3D> c{0UL, r, -hz, hz};

    ASSERT_FLOAT_EQ(c[cylinder3D::e_r], r);
    ASSERT_FLOAT_EQ(c[cylinder3D::e_n_half_z], -hz);
    ASSERT_FLOAT_EQ(c[cylinder3D::e_p_half_z], hz);

    ASSERT_TRUE(c.is_inside(p3_in) == intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside(p3_edge) == intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside(p3_out) == intersection::status::e_outside);
    ASSERT_TRUE(c.is_inside(p3_off) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside(p3_out, 0.6) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = c.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < 2; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j && i < decltype(c)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}

/// This tests the basic functionality of a 2D cylinder
TEST(mask, cylinder2D) {
    using point_t = typename mask<cylinder2D<>>::loc_point_t;

    point_t p2_in = {r, -1.};
    point_t p2_edge = {r, hz};
    point_t p2_out = {3.5, 4.5};

    // Test radius to be on surface, too
    mask<cylinder2D<>> c{0UL, r, -hz, hz};

    ASSERT_FLOAT_EQ(c[cylinder2D<>::e_r], r);
    ASSERT_FLOAT_EQ(c[cylinder2D<>::e_n_half_z], -hz);
    ASSERT_FLOAT_EQ(c[cylinder2D<>::e_p_half_z], hz);

    ASSERT_TRUE(c.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(c.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c.is_inside(p2_out, 0.6) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = c.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < 2; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j && i < decltype(c)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
}