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

/// This tests the basic functionality of a stereo annulus
TEST(mask, annulus2D) {
    using point_t = typename mask<annulus2D<>>::loc_point_t;

    constexpr scalar minR{7.2 * unit_constants::mm};
    constexpr scalar maxR{12.0 * unit_constants::mm};
    constexpr scalar minPhi{0.74195};
    constexpr scalar maxPhi{1.33970};
    point_t offset = {-2., 2.};

    // points in cartesian module frame
    point_t p2_in = {7., 7.};
    point_t p2_out1 = {5., 5.};
    point_t p2_out2 = {10., 3.};
    point_t p2_out3 = {10., 10.};
    point_t p2_out4 = {4., 10.};

    auto toStripFrame = [&](const point_t& xy) -> point_t {
        auto shifted = xy + offset;
        scalar r{getter::perp(shifted)};
        scalar phi{getter::phi(shifted)};
        return point_t{r, phi};
    };

    mask<annulus2D<>> ann2{0UL,    minR,      maxR,      minPhi,
                           maxPhi, offset[0], offset[1], 0.f};

    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_min_r], scalar{7.2});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_max_r], scalar{12.0});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_min_phi_rel], scalar{0.74195});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_max_phi_rel], scalar{1.33970});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_shift_x], scalar{-2.0});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_shift_y], scalar{2.0});
    ASSERT_FLOAT_EQ(ann2[annulus2D<>::e_average_phi], scalar{0.});

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
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out1), 1.3) ==
                intersection::status::e_inside);
    ASSERT_TRUE(ann2.is_inside(toStripFrame(p2_out4), 0.07) ==
                intersection::status::e_inside);

    // Check projection matrix
    const auto proj = ann2.projection_matrix<e_bound_size>();
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < decltype(ann2)::shape::meas_dim; j++) {
            if (i == j) {
                ASSERT_EQ(getter::element(proj, i, j), 1);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0);
            }
        }
    }

    // Test to_measurement function
    test_param param(1, 2);

    const auto meas = ann2.get_shape().to_measurement(param, {-3, 2});
    ASSERT_FLOAT_EQ(meas[0], -2.);
    ASSERT_FLOAT_EQ(meas[1], 4.);
}
