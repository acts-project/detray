/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

namespace {

constexpr scalar tol{1e-7f};

// 50 mm wire with 1 mm radial cell size
constexpr scalar cell_size{1.f * unit<scalar>::mm};
constexpr scalar hz{50.f * unit<scalar>::mm};

}  // anonymous namespace

/// This tests the basic functionality of a line with a radial cross section
TEST(mask, line_radial_cross_sect) {
    using point_t = typename mask<line<>>::loc_point_t;

    const point_t ln_in{0.09f, 0.5f};
    const point_t ln_edge{1.f, 50.f};
    const point_t ln_out1{1.2f, 0.f};
    const point_t ln_out2{0.09f, -51.f};

    const mask<line<>> ln{0u, cell_size, hz};

    ASSERT_NEAR(ln[line<>::e_cross_section], 1.f * unit<scalar>::mm, tol);
    ASSERT_NEAR(ln[line<>::e_half_z], 50.f * unit<scalar>::mm, tol);

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_out1) == intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out2) == intersection::status::e_outside);

    // Check projection matrix
    const auto proj = ln.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(ln)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }
}

/// This tests the basic functionality of a line with a square cross section
TEST(mask, line_square_cross_sect) {
    using point_t = typename mask<line<true>>::loc_point_t;

    const point_t ln_in{0.9f, 0.9f, 0.f};
    const point_t ln_edge{1.f, 1.f, 0.f};
    const point_t ln_out{1.1f, 0.f, 0.f};

    // 50 mm wire with 1 mm square cell sizes
    const mask<line<true>> ln{0u, cell_size, hz};

    ASSERT_NEAR(ln[line<>::e_cross_section], 1.f * unit<scalar>::mm, tol);
    ASSERT_NEAR(ln[line<>::e_half_z], 50.f * unit<scalar>::mm, tol);

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, 1e-5f) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, -1e-5f) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out) == intersection::status::e_outside);

    // Check projection matrix
    const auto proj = ln.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(ln)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }
}