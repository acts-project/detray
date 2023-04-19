/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/units.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/intersectors/helix_cylinder_intersector.hpp"
#include "tests/common/tools/intersectors/helix_plane_intersector.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace detray;

using transform3 = __plugin::transform3<detray::scalar>;
using vector3 = typename transform3::vector3;
using surface_handle_t = dindex;

constexpr const scalar tol = scalar(1e-4);

// dummy surface
constexpr surface_handle_t sf_handle{};

TEST(tools, helix_plane_intersector) {

    // Track defined on origin point
    const free_track_parameters<transform3> free_trk(
        {0.f, 0.f, 0.f}, 0.f, {0.1f * unit<scalar>::GeV, 0.f, 0.f}, -1.f);

    // Magnetic field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};

    const detail::helix<transform3> hlx(free_trk, &B);

    const scalar path = 10.f * unit<scalar>::cm;

    // Transform translation vector
    const auto trl = hlx(path);

    // Surface normal vector
    const auto w = hlx.dir(path);

    // z axis
    const vector3 z_axis{0.f, 0.f, 1.f};

    // Vector on the surface
    const auto v = vector::cross(z_axis, w);

    // Transform matrix
    const transform3 trf(trl, w, v);

    // Rectangle surface
    const detray::mask<detray::rectangle2D<>> rectangle{
        0u, 10.f * unit<scalar>::cm, 10.f * unit<scalar>::cm};

    const detail::helix_plane_intersector<transform3> hpi;

    // Get the intersection on the next surface
    const auto is = hpi(hlx, sf_handle, rectangle, trf);

    // Check the values
    EXPECT_NEAR(is.path, path, tol);
    EXPECT_NEAR(is.p2[0], 0.f, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_NEAR(is.p3[0], trl[0], tol);
    EXPECT_NEAR(is.p3[1], trl[1], tol);
    EXPECT_NEAR(is.p3[2], trl[2], tol);
}

TEST(tools, helix_cylinder_intersector) {

    // Track defined on origin point
    const free_track_parameters<transform3> free_trk(
        {0.f, 0.f, 0.f}, 0.f, {0.1f * unit<scalar>::GeV, 0.f, 0.f}, -1.f);

    // Magnetic field
    const vector3 B{0.f, 0.f, 1.f * unit<scalar>::T};

    const detail::helix<transform3> hlx(free_trk, &B);

    const scalar path = 10.f * unit<scalar>::cm;

    // Transform translation vector
    const auto trl = hlx(path);

    // Vector on the surface
    const auto v = hlx.dir(path);

    // z axis
    const vector3 z_axis{0.f, 0.f, 1.f};

    // Transform matrix
    const transform3 trf(trl, z_axis, v);

    const scalar c_rad = 5.f * unit<scalar>::cm;

    // Cylinder surface (5 cm radius)
    const detray::mask<detray::cylinder2D<>> cylinder{
        0u, c_rad, 10.f * unit<scalar>::cm, 10.f * unit<scalar>::cm};

    const detail::helix_cylinder_intersector<transform3> hci;

    // Get the intersection on the next surface
    const auto is = hci(hlx, sf_handle, cylinder, trf)[0];
    const auto pos = hlx.pos(is.path);
    const auto loc = pos - trl;
    const auto phi =
        std::acos(vector::dot(v, loc) / (getter::norm(v) * getter::norm(loc)));

    EXPECT_NEAR(is.p2[0], c_rad * phi, tol);
    EXPECT_NEAR(is.p2[1], 0.f, tol);
    EXPECT_NEAR(is.p3[0], pos[0], tol);
    EXPECT_NEAR(is.p3[1], pos[1], tol);
    EXPECT_NEAR(is.p3[2], pos[2], tol);
}
