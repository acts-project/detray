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

constexpr scalar tol{1e-7f};

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

    auto toStripFrame = [&](const point_t& xy) -> point_t {
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
}
