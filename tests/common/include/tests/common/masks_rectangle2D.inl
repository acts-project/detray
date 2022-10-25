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

struct test_param {
    using point2 = __plugin::point2<scalar>;

    test_param(scalar loc_0, scalar loc_1) {
        loc[0] = loc_0;
        loc[1] = loc_1;
    }

    point2 loc;
    point2 local() const { return loc; }
};

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
    for (std::size_t i = 0; i < 2; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j && i < decltype(r2)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    test_param param(1, 2);

    const auto meas = r2.get_shape().to_measurement(param, {-3, 2});
    ASSERT_FLOAT_EQ(meas[0], -2.);
    ASSERT_FLOAT_EQ(meas[1], 4.);
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
