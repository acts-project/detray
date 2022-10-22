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
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(ln)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }
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
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(ln)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    struct test_param {
        using point2 = point_t;
        point_t loc;
        point_t local() const { return loc; }
    } param_1, param_2;
    param_1.loc = {1, 2};
    param_2.loc = {2.5, 3};

    const auto meas_1 = ln.get_shape().to_measurement(param_1, {-3, 2});
    const auto meas_2 = ln.get_shape().to_measurement(param_2, {1, -4});
    ASSERT_EQ(meas_1, point_t({0, 4}));
    ASSERT_EQ(meas_2, point_t({3.5, -1}));
}