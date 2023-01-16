/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/intersection.hpp"

// Google test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <cmath>
#include <limits>

/// @note __plugin has to be defined with a preprocessor command

using namespace detray;

using transform3 = __plugin::transform3<detray::scalar>;
using point2 = __plugin::point2<scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;

/// Define mask types
enum mask_ids : unsigned int {
    e_unmasked = 0,
};

/// Define material types
enum material_ids : unsigned int {
    e_slab = 0,
};

using mask_link_t = dtyped_index<mask_ids, dindex>;
using material_link_t = dtyped_index<material_ids, dindex>;
using surface_t = surface<mask_link_t, material_link_t, transform3>;

constexpr scalar epsilon = std::numeric_limits<scalar>::epsilon();

// This tests the construction of a surface_base object
TEST(ALGEBRA_PLUGIN, surface) {
    // Preparatioon work, create a transform
    vector3 z = vector::normalize(vector3{3., 2., 1.});
    vector3 x = vector::normalize(vector3{2., -3., 0.});
    point3 t{2., 3., 4.};
    transform3 trf(t, z, x);

    mask_link_t mask_id{mask_ids::e_unmasked, 0};
    material_link_t material_id{material_ids::e_slab, 0};
    surface_t s(std::move(trf), mask_id, material_id, -1, false,
                surface_id::e_sensitive);
}

// This tests the construction of a intresection
TEST(ALGEBRA_PLUGIN, intersection) {

    using intersection_t = line_plane_intersection<surface_t, transform3>;

    intersection_t i0 = {
        intersection::status::e_outside,
        intersection::direction::e_along,
        2.f,
        1.f,
        1UL,
        surface_t{},
        point3{0.3f, 0.5f, 0.7f},
        point2{0.2f, 0.4f},
    };

    intersection_t i1 = {
        intersection::status::e_inside,
        intersection::direction::e_opposite,
        1.7f,
        -1.f,
        0UL,
        surface_t{},
        point3{0.2f, 0.3f, 0.f},
        point2{0.2f, 0.4f},
    };

    intersection_t invalid;
    ASSERT_TRUE(invalid.status == intersection::status::e_undefined);

    dvector<intersection_t> intersections = {invalid, i0, i1};
    std::sort(intersections.begin(), intersections.end());

    ASSERT_NEAR(intersections[0].path, 1.7f, epsilon);
    ASSERT_NEAR(intersections[1].path, 2.f, epsilon);
    ASSERT_TRUE(std::isinf(intersections[2].path));
}
