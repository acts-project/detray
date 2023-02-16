/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;
using namespace __plugin;

constexpr scalar tol{1e-5f};

/// This tests the basic functionality of a stereo annulus
TEST(mask, annulus2D) {
    using point_t = typename mask<annulus2D<>>::loc_point_t;

    constexpr scalar minR{7.2f * unit<scalar>::mm};
    constexpr scalar maxR{12.0f * unit<scalar>::mm};
    constexpr scalar minPhi{0.74195f};
    constexpr scalar maxPhi{1.33970f};
    point_t offset = {-2.f, 2.f};

    // points in cartesian module frame
    point_t p2_in = {7.f, 7.f};
    point_t p2_out1 = {5.f, 5.f};
    point_t p2_out2 = {10.f, 3.f};
    point_t p2_out3 = {10.f, 10.f};
    point_t p2_out4 = {4.f, 10.f};

    auto toStripFrame = [&offset](const point_t& xy) -> point_t {
        auto shifted = xy + offset;
        scalar r{getter::perp(shifted)};
        scalar phi{getter::phi(shifted)};
        return point_t{r, phi};
    };

    mask<annulus2D<>> ann2{0u,     minR,      maxR,      minPhi,
                           maxPhi, offset[0], offset[1], 0.f};

    ASSERT_NEAR(ann2[annulus2D<>::e_min_r], 7.2f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_max_r], 12.0f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_min_phi_rel], 0.74195f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_max_phi_rel], 1.33970f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_shift_x], -2.0f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_shift_y], 2.0f, tol);
    ASSERT_NEAR(ann2[annulus2D<>::e_average_phi], 0.f, tol);

    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_in)) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out2)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out3)) ==
                intersection::status::e_outside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4)) ==
                intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1), 1.3f) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4), 0.07f) ==
                intersection::status::e_inside);

    // Check projection matrix
    const auto proj = ann2.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(ann2)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check corner positions
    std::array<scalar, 8> c = ann2.get_shape().corners(ann2.values());
    for (unsigned int i{0u}; i < 8u; i += 2u) {
        // Transform to local cartesian beam system
        const scalar loc_x{c[i] * math_ns::cos(c[i + 1]) -
                           ann2.values()[annulus2D<>::e_shift_x]};
        const scalar loc_y{c[i] * math_ns::sin(c[i + 1]) -
                           ann2.values()[annulus2D<>::e_shift_y]};

        // Inner points
        if (i < 4u) {
            EXPECT_NEAR(std::hypot(loc_x, loc_y), minR, tol)
                << "point " << i << ": loc_x: " << loc_x << ", loc_y: " << loc_y
                << std::endl;
        }
        // Outer points
        else {
            EXPECT_NEAR(std::hypot(loc_x, loc_y), maxR, tol)
                << "point " << i << ": loc_x: " << loc_x << ", loc_y: " << loc_y
                << std::endl;
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = ann2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], 3.8954f - envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], 2.39186f - envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], 10.50652f + envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], 10.89317f + envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}
