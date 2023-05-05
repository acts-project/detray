/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/intersection.hpp"
#include "detray/test/types.hpp"
#include "detray/utils/invalid_values.hpp"

// Google test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <cmath>
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

using transform3 = test::transform3;
using point2 = test::point2;
using vector3 = test::vector3;
using point3 = test::point3;

using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;
using surface_t = surface<mask_link_t, material_link_t, transform3>;

/// Test the typed index
TEST(core, typed_index) {

    using index_t = dtyped_index<mask_ids, unsigned int>;
    auto ti = index_t{};

    // Check a empty barcode
    EXPECT_EQ(ti.id(), static_cast<unsigned int>((1u << 4) - 1u));
    EXPECT_EQ(ti.index(), static_cast<unsigned int>((1u << 28) - 1u));

    ti.set_id(mask_ids::e_unmasked).set_index(42u);

    // Check the values after setting them
    EXPECT_EQ(ti.id(), mask_ids::e_unmasked);
    EXPECT_EQ(ti.index(), 42u);

    // Check invalid link
    EXPECT_FALSE(ti.is_invalid());
    ti.set_id(static_cast<index_t::id_type>((1u << 4) - 1u));
    EXPECT_TRUE(ti.is_invalid());
    ti.set_id(mask_ids::e_unmasked);
    EXPECT_FALSE(ti.is_invalid());
    ti.set_index((1u << 30) - 1u);
    EXPECT_TRUE(ti.is_invalid());
}

// This tests the construction of a surface descriptor object
TEST(core, surface) {

    mask_link_t mask_id{mask_ids::e_unmasked, 0u};
    material_link_t material_id{material_ids::e_slab, 0u};

    surface<mask_link_t, material_link_t> s(
        dindex_invalid, mask_id, material_id, dindex_invalid, dindex_invalid,
        surface_id::e_sensitive);
}

// This tests the construction of a intresection
GTEST_TEST(detray_core, intersection) {

    using intersection_t = intersection2D<surface_t, transform3>;

    intersection_t i0 = {surface_t{},
                         point3{0.2f, 0.4f, 0.f},
                         2.f,
                         1.f,
                         1u,
                         intersection::status::e_outside,
                         intersection::direction::e_along};

    intersection_t i1 = {
        surface_t{},
        point3{0.2f, 0.4f, 0.f},
        1.7f,
        -1.f,
        0u,
        intersection::status::e_inside,
        intersection::direction::e_opposite,
    };

    intersection_t invalid;
    ASSERT_TRUE(invalid.status == intersection::status::e_undefined);

    dvector<intersection_t> intersections = {invalid, i0, i1};
    std::sort(intersections.begin(), intersections.end());

    ASSERT_NEAR(intersections[0].path, 1.7f, tol);
    ASSERT_NEAR(intersections[1].path, 2.f, tol);
    ASSERT_TRUE(is_invalid_value(intersections[2].path));
}
