/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <cassert>
#include <type_traits>

using namespace detray;
using point3_t = __plugin::point3<detray::scalar>;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of an unbounded rectangle shape
TEST(mask, unbounded) {
    using transform3_t = __plugin::transform3<scalar>;

    using shape_t = rectangle2D<>;
    using unbounded_t = unbounded<shape_t>;

    constexpr scalar h{20.f * unit<scalar>::mm};

    mask<unbounded_t> u{0u, h, h};

    // Test local typedefs
    static_assert(std::is_same_v<unbounded_t::shape, shape_t>,
                  "incorrect shape");
    static_assert(std::is_same_v<unbounded_t::boundaries, shape_t::boundaries>,
                  "incorrect boundaries");
    static_assert(
        std::is_same_v<unbounded_t::template local_frame_type<transform3_t>,
                       shape_t::template local_frame_type<transform3_t>>,
        "incorrect local frame");
    static_assert(
        std::is_same_v<unbounded_t::template intersector_type<transform3_t>,
                       shape_t::template intersector_type<transform3_t>>,
        "incorrect intersector");

    // Test static members
    EXPECT_TRUE(unbounded_t::name == "unbounded rectangle2D");
    EXPECT_TRUE(unbounded_t::meas_dim == 2u);

    // Test boundary check
    typename mask<unbounded_t>::point3_t p2 = {0.5f, -9.f, 0.f};
    ASSERT_TRUE(u.is_inside(p2, 0.f) == intersection::status::e_inside);

    // Check projection matrix
    const auto proj = u.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
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
    const auto loc_bounds = u.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(h + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(h + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (h + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (h + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}
