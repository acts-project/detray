/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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

namespace {

// 50 mm wire with 1 mm radial cell size
constexpr scalar cell_size{1. * unit_constants::mm};
constexpr scalar hz{50. * unit_constants::mm};

}  // anonymous namespace

/// This tests the basic functionality of a line with a radial cross section
TEST(mask, line_radial_cross_sect) {
    using point_t = typename mask<line<>>::loc_point_t;

    const point_t ln_in{0.09, 0.5};
    const point_t ln_edge{1., 50.};
    const point_t ln_out1{1.2, 0};
    const point_t ln_out2{0.09, -51.};

    const mask<line<>> ln{0UL, cell_size, hz};

    ASSERT_FLOAT_EQ(ln[line<>::e_cross_section],
                    scalar{1. * unit_constants::mm});
    ASSERT_FLOAT_EQ(ln[line<>::e_half_z], scalar{50. * unit_constants::mm});

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_out1) == intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out2) == intersection::status::e_outside);

    // Check projection matrix
    const auto proj = ln.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < 2; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j && i < decltype(ln)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    test_param param_1(1, 2);
    test_param param_2(2.5, 3);

    const auto meas_1 = ln.get_shape().to_measurement(param_1, {-3, 2});
    const auto meas_2 = ln.get_shape().to_measurement(param_2, {1, -4});
    ASSERT_FLOAT_EQ(meas_1[0], 0.);
    ASSERT_FLOAT_EQ(meas_1[1], 4.);
    ASSERT_FLOAT_EQ(meas_2[0], 3.5);
    ASSERT_FLOAT_EQ(meas_2[1], -1);
}

/// This tests the basic functionality of a line with a square cross section
TEST(mask, line_square_cross_sect) {
    using point_t = typename mask<line<true>>::loc_point_t;

    const point_t ln_in{0.9, 0.9, 0};
    const point_t ln_edge{1., 1., 0.};
    const point_t ln_out{1.1, 0., 0};

    // 50 mm wire with 1 mm square cell sizes
    const mask<line<true>> ln{0UL, cell_size, hz};

    ASSERT_FLOAT_EQ(ln[line<>::e_cross_section],
                    scalar{1. * unit_constants::mm});
    ASSERT_FLOAT_EQ(ln[line<>::e_half_z], scalar{50. * unit_constants::mm});

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, 1e-5) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, -1e-5) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out) == intersection::status::e_outside);

    // Check projection matrix
    const auto proj = ln.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < 2; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            if (i == j && i < decltype(ln)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    test_param param_1(1, 2);
    test_param param_2(2.5, 3);

    const auto meas_1 = ln.get_shape().to_measurement(param_1, {-3, 2});
    const auto meas_2 = ln.get_shape().to_measurement(param_2, {1, -4});
    ASSERT_FLOAT_EQ(meas_1[0], 0.);
    ASSERT_FLOAT_EQ(meas_1[1], 4.);
    ASSERT_FLOAT_EQ(meas_2[0], 3.5);
    ASSERT_FLOAT_EQ(meas_2[1], -1);
}