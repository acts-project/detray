/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/masks/masks.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

using namespace detray;
using point3_t = test::point3;
using transform3_t = test::transform3;

/// This tests the basic functionality of an unmasked plane
GTEST_TEST(detray_masks, unmasked) {
    point3_t p2 = {0.5f, -9.f, 0.f};

    mask<unmasked> u{};

    ASSERT_TRUE(u.is_inside(p2, 0.f) == intersection::status::e_inside);

    // Dummy bound track parameter
    bound_track_parameters<transform3_t> bound_params;

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = u.local_min_bounds(envelope);
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_min_x]));
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_min_y]));
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_min_z]));
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_max_x]));
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_max_y]));
    ASSERT_TRUE(std::isinf(loc_bounds[cuboid3D<>::e_max_z]));
}
