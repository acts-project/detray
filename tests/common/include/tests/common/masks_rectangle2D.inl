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

constexpr scalar hx{1. * unit_constants::mm};
constexpr scalar hy{9.3 * unit_constants::mm};
constexpr scalar hz{0.5 * unit_constants::mm};

/// This tests the basic functionality of a rectangle
TEST(mask, rectangle2D) {
    using point_t = typename mask<rectangle2D<>>::loc_point_t;

    point_t p2_in = {0.5, -9.};
    point_t p2_edge = {1., 9.3};
    point_t p2_out = {1.5, -9.};

    mask<rectangle2D<>> r2{0UL, hx, hy};

    ASSERT_FLOAT_EQ(r2[rectangle2D<>::e_half_x], hx);
    ASSERT_FLOAT_EQ(r2[rectangle2D<>::e_half_y], hy);

    ASSERT_TRUE(r2.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_out, 1.) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = r2.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(r2)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    struct test_param {
        using point2 = point_t;
        point_t loc;
        point_t local() const { return loc; }
    } param;
    param.loc = {1, 2};

    const auto meas = r2.get_shape().to_measurement(param, {-3, 2});
    ASSERT_EQ(meas, point_t({-2, 4}));
}

/// This tests the basic functionality of a cuboid3D
TEST(mask, cuboid3D) {
    using point_t = typename mask<cuboid3D>::loc_point_t;

    point_t p2_in = {0.5, 8.0, -0.4};
    point_t p2_edge = {1., 9.3, 0.5};
    point_t p2_out = {1.5, -9., 0.55};

    mask<cuboid3D> c3{0UL, hx, hy, hz};

    ASSERT_FLOAT_EQ(c3[cuboid3D::e_half_x], hx);
    ASSERT_FLOAT_EQ(c3[cuboid3D::e_half_y], hy);
    ASSERT_FLOAT_EQ(c3[cuboid3D::e_half_z], hz);

    ASSERT_TRUE(c3.is_inside(p2_in) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_edge) == intersection::status::e_inside);
    ASSERT_TRUE(c3.is_inside(p2_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(c3.is_inside(p2_out, 1.) == intersection::status::e_inside);
}
