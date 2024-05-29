/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/detail/indexing.hpp"
#include "detray/geometry/detail/surface_descriptor.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/test/common/types.hpp"
#include "detray/utils/invalid_values.hpp"

// Google test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

namespace {

/// Define mask types
enum mask_ids : unsigned int {
    e_unmasked = 0u,
};

/// Define material types
enum material_ids : unsigned int {
    e_slab = 0u,
};

constexpr scalar tol{std::numeric_limits<scalar>::epsilon()};

}  // namespace

using algebra_t = test::algebra;
using point2 = test::point2;
using vector3 = test::vector3;
using point3 = test::point3;

using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;
using surface_t = surface_descriptor<mask_link_t, material_link_t, algebra_t>;

// This tests the construction of a intresection
GTEST_TEST(detray_intersection, intersection2D) {

    using intersection_t = intersection2D<surface_t, algebra_t>;

    intersection_t i0 = {
        surface_t{}, point3{0.2f, 0.4f, 0.f}, 2.f, 1.f, 1u, false, true};

    intersection_t i1 = {
        surface_t{}, point3{0.2f, 0.4f, 0.f}, 1.7f, -1.f, 0u, true, false,
    };

    intersection_t invalid{};
    ASSERT_FALSE(invalid.status);

    dvector<intersection_t> intersections = {invalid, i0, i1};
    std::sort(intersections.begin(), intersections.end());

    ASSERT_NEAR(intersections[0].path, 1.7f, tol);
    ASSERT_NEAR(intersections[1].path, 2.f, tol);
    ASSERT_TRUE(detail::is_invalid_value(intersections[2].path));
}
