/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unmasked.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"
#include "tests/common/tools/intersectors/helix_line_intersector.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s)
#include <cmath>
#include <limits>

using namespace detray;

namespace {

// Three-dimensional definitions
using transform3_t = __plugin::transform3<detray::scalar>;
using vector3 = __plugin::vector3<detray::scalar>;
using point3 = __plugin::point3<detray::scalar>;
using helix_t = detray::detail::helix<transform3_t>;
using intersection_t = intersection2D_point<surface<>, transform3_t>;

constexpr auto not_defined{detail::invalid_value<scalar>()};
constexpr scalar tol{1e-4f};

// z axis
const vector3 z_axis{0.f, 0.f, 1.f};

// Track defined on origin point
const free_track_parameters<transform3_t> free_trk(
    {0.f, 0.f, 0.f}, 0.f, {0.1f * unit<scalar>::GeV, 0.f, 0.f}, -1.f);

// Magnetic field
const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};
const vector3 B_0{0.f * unit<scalar>::T, 0.f * unit<scalar>::T,
                  tol* unit<scalar>::T};

// Test helix
const helix_t hlx(free_trk, &B);

// Path along the helix
const scalar path = 10.f * unit<scalar>::cm;

// Transform translation vector
const vector3 trl = hlx(path);

// Surface normal vector
const vector3 w = hlx.dir(path);

}  // anonymous namespace

/// This defines the local frame test suite
TEST(tools, helix_plane_intersector_no_bfield) {
    // Create a shifted plane
    const transform3_t shifted(vector3{3.f, 2.f, 10.f});

    // Test helix
    const point3 pos{2.f, 1.f, 0.f};
    const vector3 mom{0.f, 0.f, 1.f};
    const detail::helix<transform3_t> h({pos, 0.f, mom, -1.f}, &B_0);

    // The same test but bound to local frame
    detail::helix_plane_intersector<intersection_t> pi;
    mask<unmasked> unmasked_bound{};
    const auto hit_bound = pi(h, surface<>{}, unmasked_bound, shifted);

    ASSERT_TRUE(hit_bound.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound.p3[0], 2.f, tol);
    ASSERT_NEAR(hit_bound.p3[1], 1.f, tol);
    ASSERT_NEAR(hit_bound.p3[2], 10.f, tol);
    // Local intersection information
    ASSERT_NEAR(hit_bound.p2[0], -1.f, tol);
    ASSERT_NEAR(hit_bound.p2[1], -1.f, tol);
    // Incidence angle
    ASSERT_TRUE(is_invalid_value(hit_bound.cos_incidence_angle));

    // The same test but bound to local frame & masked - inside
    mask<rectangle2D<>> rect_for_inside{0u, 3.f, 3.f};
    const auto hit_bound_inside = pi(h, surface<>{}, rect_for_inside, shifted);
    ASSERT_TRUE(hit_bound_inside.status == intersection::status::e_inside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_inside.p3[0], 2.f, tol);
    ASSERT_NEAR(hit_bound_inside.p3[1], 1.f, tol);
    ASSERT_NEAR(hit_bound_inside.p3[2], 10.f, tol);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_inside.p2[0], -1.f, tol);
    ASSERT_NEAR(hit_bound_inside.p2[1], -1.f, tol);

    // The same test but bound to local frame & masked - outside
    mask<rectangle2D<>> rect_for_outside{0u, 0.5f, 3.5f};
    const auto hit_bound_outside =
        pi(h, surface<>{}, rect_for_outside, shifted);
    ASSERT_TRUE(hit_bound_outside.status == intersection::status::e_outside);
    // Global intersection information - unchanged
    ASSERT_NEAR(hit_bound_outside.p3[0], 2.f, tol);
    ASSERT_NEAR(hit_bound_outside.p3[1], 1.f, tol);
    ASSERT_NEAR(hit_bound_outside.p3[2], 10.f, tol);
    // Local intersection infoimation - unchanged
    ASSERT_NEAR(hit_bound_outside.p2[0], -1.f, tol);
    ASSERT_NEAR(hit_bound_outside.p2[1], -1.f, tol);
}

/// Test the intersection between a helical trajectory and a plane
TEST(tools, helix_plane_intersector) {

    // Vector on the surface
    const vector3 v = vector::cross(z_axis, w);

    // Transform matrix
    const transform3_t trf(trl, w, v);

    // Rectangle surface
    const mask<rectangle2D<>> rectangle{0u, 10.f * unit<scalar>::cm,
                                        10.f * unit<scalar>::cm};

    const detail::helix_plane_intersector<intersection_t> hpi;

    // Get the intersection on the next surface
    const auto is = hpi(hlx, surface<>{}, rectangle, trf, tol);

    // Check the values
    EXPECT_TRUE(is.status == intersection::status::e_inside);
    EXPECT_NEAR(is.path, path, tol);
    EXPECT_NEAR(is.p2[0], 0.f, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_NEAR(is.p3[0], trl[0], tol);
    EXPECT_NEAR(is.p3[1], trl[1], tol);
    EXPECT_NEAR(is.p3[2], trl[2], tol);
}

/// This checks the closest solution of a helix-cylinder intersection
TEST(tools, helix_cylinder_intersector_no_bfield) {

    const scalar r{4.f * unit<scalar>::mm};
    const scalar hz{10.f * unit<scalar>::mm};

    // Create a translated cylinder and test untersection
    const transform3_t shifted(vector3{3.f, 2.f, 10.f});
    detail::helix_cylinder_intersector<intersection_t> hi;

    // Test helix
    const point3 pos{3.f, 2.f, 5.f};
    const vector3 mom{1.f, 0.f, 0.f};
    const helix_t h({pos, 0.f * unit<scalar>::s, mom, -1 * unit<scalar>::e},
                    &B_0);

    // Intersect
    mask<cylinder2D<false, detail::helix_cylinder_intersector>,
         std::uint_least16_t, transform3_t>
        cylinder{0u, r, -hz, hz};
    const auto hits_bound = hi(h, surface<>{}, cylinder, shifted, tol);

    // No magnetic field, so the solutions must be the same as for a ray

    // second intersection lies in front of the track
    EXPECT_TRUE(hits_bound[0].status == intersection::status::e_inside);
    EXPECT_TRUE(hits_bound[0].direction == intersection::direction::e_opposite);
    EXPECT_NEAR(hits_bound[0].p3[0], -1.f, tol);
    EXPECT_NEAR(hits_bound[0].p3[1], 2.f, tol);
    EXPECT_NEAR(hits_bound[0].p3[2], 5.f, tol);
    ASSERT_TRUE(hits_bound[0].p2[0] != not_defined &&
                hits_bound[0].p2[1] != not_defined);
    // p2[0] = r * phi : 180deg in the opposite direction with r = 4
    EXPECT_NEAR(hits_bound[0].p2[0], 4.f * M_PI, tol);
    EXPECT_NEAR(hits_bound[0].p2[1], -5.f, tol);
    EXPECT_TRUE(is_invalid_value(hits_bound[0].cos_incidence_angle));

    // first intersection lies behind the track
    EXPECT_TRUE(hits_bound[1].status == intersection::status::e_inside);
    EXPECT_TRUE(hits_bound[1].direction == intersection::direction::e_along);
    EXPECT_NEAR(hits_bound[1].p3[0], 7.f, tol);
    EXPECT_NEAR(hits_bound[1].p3[1], 2.f, tol);
    EXPECT_NEAR(hits_bound[1].p3[2], 5.f, tol);
    ASSERT_TRUE(hits_bound[1].p2[0] != not_defined &&
                hits_bound[1].p2[1] != not_defined);
    EXPECT_NEAR(hits_bound[1].p2[0], 0.f, tol);
    EXPECT_NEAR(hits_bound[1].p2[1], -5., tol);
    EXPECT_TRUE(is_invalid_value(hits_bound[1].cos_incidence_angle));
}

/// Test the intersection between a helical trajectory and a cylinder
TEST(tools, helix_cylinder_intersector) {

    // Transform matrix
    const transform3_t trf(trl, z_axis, w);

    // Cylinder surface (5 cm radius)
    const scalar r{4.f * unit<scalar>::cm};
    const scalar hz{10.f * unit<scalar>::cm};
    const mask<cylinder2D<>> cylinder{0u, r, -hz, hz};

    const detail::helix_cylinder_intersector<intersection_t> hci;

    // Get the intersection on the next surface
    const auto is = hci(hlx, surface<>{}, cylinder, trf, tol);

    // First solution
    const vector3 pos_near = hlx.pos(is[0].path);
    const vector3 loc_near = pos_near - trl;
    const scalar phi_near = std::acos(
        vector::dot(w, loc_near) / (getter::norm(w) * getter::norm(loc_near)));

    EXPECT_TRUE(is[0].status == intersection::status::e_inside);
    // Not precise due to helix curvature
    EXPECT_NEAR(is[0].path, path - r, 5000.f * tol);
    EXPECT_NEAR(is[0].p2[0], r * phi_near, tol);
    EXPECT_NEAR(is[0].p2[1], 0.f, tol);
    EXPECT_NEAR(is[0].p3[0], pos_near[0], tol);
    EXPECT_NEAR(is[0].p3[1], pos_near[1], tol);
    EXPECT_NEAR(is[0].p3[2], pos_near[2], tol);

    // Second solution
    const vector3 pos_far = hlx.pos(is[1].path);
    const vector3 loc_far = pos_far - trl;
    const scalar phi_far = std::acos(vector::dot(w, loc_far) /
                                     (getter::norm(w) * getter::norm(loc_far)));

    EXPECT_TRUE(is[1].status == intersection::status::e_inside);
    // Not precise due to helix curvature
    EXPECT_NEAR(is[1].path, path + r, 5000.f * tol);
    EXPECT_NEAR(is[1].p2[0], r * phi_far, tol);
    EXPECT_NEAR(is[1].p2[1], 0.f, tol);
    EXPECT_NEAR(is[1].p3[0], pos_far[0], tol);
    EXPECT_NEAR(is[1].p3[1], pos_far[1], tol);
    EXPECT_NEAR(is[1].p3[2], pos_far[2], tol);
}

/// Test the intersection between a helical trajectory and a line
TEST(tools, helix_line_intersector) {

    // Intersector object
    const detail::helix_line_intersector<intersection_t> hli;

    // Get radius of track
    const scalar R{hlx.radius()};

    // Path length for pi/4
    const scalar s0 = constant<scalar>::pi_2 * R;

    // Wire properties
    const scalar scope = 2.f * unit<scalar>::cm;
    const scalar half_z = std::numeric_limits<scalar>::max();

    // Straw wire
    const mask<line<false>> straw_wire{0u, scope, half_z};

    // Cell wire
    const mask<line<true>> cell_wire{0u, scope, half_z};

    // Offset to shift the translation of transform matrix
    const scalar offset = 1.f * unit<scalar>::cm;

    //---------------------
    // Forward direction
    //---------------------

    // Reference point for transform matrix
    const point3 r0_fw = hlx.pos(s0);

    // Translation is shifted from reference point
    const point3 trl_fw = r0_fw + vector3{offset, 0.f, 0.f};

    // Transform matrix
    const transform3_t trf_fw(trl_fw, z_axis, hlx.dir(s0));

    // Get the intersection on the next surface
    auto is = hli(hlx, surface<>{}, straw_wire, trf_fw, tol);

    EXPECT_NEAR(is.path, s0, tol);
    // track (helix) is at the left side w.r.t wire
    EXPECT_NEAR(is.p2[0], offset, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_EQ(is.direction, intersection::direction::e_along);

    // Get the intersection on the next surface
    is = hli(hlx, surface<>{}, cell_wire, trf_fw, tol);

    EXPECT_NEAR(is.path, s0, tol);
    // track (helix) is at the left side w.r.t wire
    EXPECT_NEAR(is.p2[0], offset, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_EQ(is.direction, intersection::direction::e_along);

    //---------------------
    // Backward direction
    //---------------------

    // Reference point for transform matrix
    const point3 r0_bw = hlx.pos(-s0);

    // Translation is shifted from reference point
    const point3 trl_bw = r0_bw + vector3{offset, 0.f, 0.f};

    // Transform matrix
    const transform3_t trf_bw(trl_bw, z_axis, hlx.dir(-s0));

    // Get the intersection on the next surface
    is = hli(hlx, surface<>{}, straw_wire, trf_bw, tol);

    EXPECT_NEAR(is.path, -s0, tol);
    // track (helix) is at the right side w.r.t wire
    EXPECT_NEAR(is.p2[0], -offset, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_EQ(is.direction, intersection::direction::e_opposite);

    // Get the intersection on the next surface
    is = hli(hlx, surface<>{}, cell_wire, trf_bw, tol);

    EXPECT_NEAR(is.path, -s0, tol);
    // track (helix) is at the right side w.r.t wire
    EXPECT_NEAR(is.p2[0], -offset, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_EQ(is.status, intersection::status::e_inside);
    EXPECT_EQ(is.direction, intersection::direction::e_opposite);
}