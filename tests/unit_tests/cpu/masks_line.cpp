/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

using namespace detray;
using point3_t = test::point3;
using transform3_t = test::transform3;

namespace {

constexpr scalar tol{1e-7f};

// 50 mm wire with 1 mm radial cell size
constexpr scalar cell_size{1.f * unit<scalar>::mm};
constexpr scalar hz{50.f * unit<scalar>::mm};

}  // anonymous namespace

/// This tests the basic functionality of a line with a radial cross section
GTEST_TEST(detray_masks, line_radial_cross_sect) {
    using point_t = point3_t;

    const point_t ln_in{0.09f, 0.5f, 0.f};
    const point_t ln_edge{1.f, 50.f, 0.f};
    const point_t ln_out1{1.2f, 0.f, 0.f};
    const point_t ln_out2{0.09f, -51.f, 0.f};

    const mask<line<>> ln{0u, cell_size, hz};

    ASSERT_NEAR(ln[line<>::e_cross_section], 1.f * unit<scalar>::mm, tol);
    ASSERT_NEAR(ln[line<>::e_half_z], 50.f * unit<scalar>::mm, tol);

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_out1) == intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out2) == intersection::status::e_outside);

    // Dummy bound track parameter
    bound_track_parameters<transform3_t> bound_params;
    auto& bound_vec = bound_params.vector();
    getter::element(bound_vec, e_bound_loc0, 0u) = 1.f;

    // Check projection matrix
    const auto proj = ln.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < decltype(ln)::shape::meas_dim; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = ln.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}

/// This tests the basic functionality of a line with a square cross section
GTEST_TEST(detray_masks, line_square_cross_sect) {
    using point_t = point3_t;

    const point_t ln_in{1.1f, 0.9f, constant<scalar>::pi_4};
    const point_t ln_edge{1.f, 1.f, 0.f};
    const point_t ln_out1{1.1f, 0.f, 0.f};
    const point_t ln_out2{0.09f, -51.f, 0.f};

    // 50 mm wire with 1 mm square cell sizes
    const mask<line<true>> ln{0u, cell_size, hz};

    ASSERT_NEAR(ln[line<>::e_cross_section], 1.f * unit<scalar>::mm, tol);
    ASSERT_NEAR(ln[line<>::e_half_z], 50.f * unit<scalar>::mm, tol);

    ASSERT_TRUE(ln.is_inside(ln_in) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, 1e-5f) == intersection::status::e_inside);
    ASSERT_TRUE(ln.is_inside(ln_edge, -1e-5f) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out1) == intersection::status::e_outside);
    ASSERT_TRUE(ln.is_inside(ln_out2) == intersection::status::e_outside);

    // Dummy bound track parameter
    bound_track_parameters<transform3_t> bound_params;
    auto& bound_vec = bound_params.vector();
    getter::element(bound_vec, e_bound_loc0, 0u) = -1.f;

    // Check projection matrix
    const auto proj = ln.projection_matrix(bound_params);
    for (unsigned int i = 0u; i < decltype(ln)::shape::meas_dim; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j) {
                if (i == 1) {
                    ASSERT_EQ(getter::element(proj, i, j), -1.f);
                } else if (i == 2) {
                    ASSERT_EQ(getter::element(proj, i, j), 1.f);
                }
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = ln.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -(hz + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (cell_size + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], (hz + envelope), tol);
}