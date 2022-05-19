/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/mask_store.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/enumerate.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/test_trajectories.hpp"
#include "tests/common/tools/track_generators.hpp"
/// @note __plugin has to be defined with a preprocessor command
using namespace detray;
using namespace __plugin;

constexpr const float epsilon = 1e-5;

// This tests the base functionality of the Helix Gun
TEST(tools, helix_trajectory) {

    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    point3 pos{0., 0., 0.};
    scalar time = 0.;
    vector3 mom{1., 0., 1.};
    scalar q = -1.;

    // vertex
    free_track_parameters vertex(pos, time, mom, q);

    // magnetic field
    vector3 B{0, 0, 1 * unit_constants::T};

    scalar p_mag = getter::norm(mom);
    scalar B_mag = getter::norm(B);
    scalar pz_along = vector::dot(mom, vector::normalize(B));
    scalar pt = std::sqrt(std::pow(p_mag, 2) - std::pow(pz_along, 2));

    // helix trajectory
    helix helix_traj(vertex, &B);

    // radius of helix
    scalar R = helix_traj.radius();
    EXPECT_FLOAT_EQ(R, pt / B_mag);

    // After half turn
    point3 half_turn = helix_traj(p_mag / B_mag * M_PI);

    EXPECT_NEAR(half_turn[0], 0., R * epsilon);
    EXPECT_NEAR(half_turn[1], 2 * R, R * epsilon);
    EXPECT_NEAR(half_turn[2], pz_along / B_mag * M_PI, R * epsilon);

    // After one full turn
    point3 full_turn = helix_traj(2 * p_mag / B_mag * M_PI);

    EXPECT_NEAR(full_turn[0], 0., R * epsilon);
    EXPECT_NEAR(full_turn[1], 0., R * epsilon);
    EXPECT_NEAR(full_turn[2], 2 * pz_along / B_mag * M_PI, R * epsilon);
}

TEST(tools, helix_trajectory_small_pT) {

    using vector3 = vector3<scalar>;
    using point3 = point3<scalar>;

    point3 pos{0., 0., 0.};
    scalar time = 0.;
    vector3 mom{0., 0., 1.};
    scalar q = -1.;

    // vertex
    free_track_parameters vertex(pos, time, mom, q);

    // magnetic field
    vector3 B{0, 0, 1 * unit_constants::T};

    // helix trajectory
    helix helix_traj(vertex, &B);

    // After 10 mm
    scalar path_length = 10;
    point3 helix_pos = helix_traj(path_length);
    point3 true_pos = pos + path_length * vector::normalize(mom);

    EXPECT_FLOAT_EQ(true_pos[0], helix_pos[0]);
    EXPECT_FLOAT_EQ(true_pos[1], helix_pos[1]);
    EXPECT_FLOAT_EQ(true_pos[2], helix_pos[2]);
}

/// Re-use the intersection kernel test for particle gun
TEST(tools, helix_gun) {
    using vector3 = __plugin::vector3<scalar>;
    using point3 = __plugin::point3<scalar>;
    using transform3 = __plugin::transform3<scalar>;

    vecmem::host_memory_resource host_mr;

    enum mask_ids : unsigned int {
        e_rectangle2 = 0,
        e_trapezoid2 = 1,
        e_annulus2 = 2,
    };

    /// Surface components:
    using edge_t = dindex;
    using source_link_t = dindex;
    /// - masks, with mask identifiers 0,1,2
    using rectangle_t =
        rectangle2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_t>;
    using trapezoid_t =
        trapezoid2<planar_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_t>;
    using annulus_t = annulus2<planar_intersector,
                               __plugin::cartesian2<detray::scalar>, edge_t>;

    using mask_defs =
        mask_registry<mask_ids, rectangle_t, trapezoid_t, annulus_t>;
    using mask_container_t = typename mask_defs::container_type<>;

    /// The Surface definition:
    /// <transform_link, volume_link, source_link, link_type_in_mask>
    using surface_t = surface<mask_defs, dindex, dindex, source_link_t>;
    using surface_container_t = dvector<surface_t>;

    // The transforms & their store
    transform3 rectangle_transform(point3{0., 0., 10.});
    transform3 trapezoid_transform(point3{0., 0., 20.});
    transform3 annulus_transform(point3{0., -20., 30.});
    static_transform_store<>::context static_context{};
    static_transform_store transform_store;
    transform_store.push_back(static_context, rectangle_transform);
    transform_store.push_back(static_context, trapezoid_transform);
    transform_store.push_back(static_context, annulus_transform);
    // The masks & their store
    mask_container_t mask_store(host_mr);
    mask_store.template add_mask<e_rectangle2>(10., 10., 0);
    mask_store.template add_mask<e_trapezoid2>(10., 20., 30., 0);
    mask_store.template add_mask<e_annulus2>(15., 55., 0.75, 1.95, 2., -2., 0.,
                                             0);
    // The surfaces and their store
    surface_t rectangle_surface(0u, {e_rectangle2, 0}, 0, 0, false);
    surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, 0, 1, false);
    surface_t annulus_surface(2u, {e_annulus2, 0}, 0, 2, false);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};

    point3 pos{0., 0., 0.};
    vector3 mom{0.01, 0.01, 10.};
    ray r(pos, 0, mom, -1);
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    0.001 * unit_constants::T};
    free_track_parameters vertex(pos, 0, mom, -1);
    helix h(vertex, &B);
    // const auto dir = h.pos(0.);

    // std::cout << "dir: " << dir[0] << ", " << dir[1] << ", " << dir[2] <<
    // std::endl;

    // Validation data
    point3 expected_rectangle{0.01, 0.01, 10.};
    point3 expected_trapezoid{0.02, 0.02, 20.};
    point3 expected_annulus{0.03, 0.03, 30.};

    std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersection - with automated dispatching via the kernel
    unsigned int it = 0;
    for (const auto &_surface : surfaces) {
        auto sfi_helix =
            particle_gun::intersect(h, _surface, transform_store, mask_store);

        ASSERT_NEAR(sfi_helix.p3[0], expected_points[it][0], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[1], expected_points[it][1], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[2], expected_points[it][2], 1e-7);
        // std::cout << sfi_helix.p3[0] << ", " << sfi_helix.p3[1] << ", " <<
        // sfi_helix.p3[2] << std::endl;

        auto sfi_ray =
            particle_gun::intersect(r, _surface, transform_store, mask_store);
        ASSERT_NEAR(sfi_ray.p3[0], expected_points[it][0], 1e-7);
        ASSERT_NEAR(sfi_ray.p3[1], expected_points[it][1], 1e-7);
        ASSERT_NEAR(sfi_ray.p3[2], expected_points[it][2], 1e-7);
        ++it;
    }
}

TEST(tools, uniform_track_generator) {

    using vector3 = __plugin::vector3<scalar>;

    constexpr std::size_t phi_steps = 5;
    constexpr std::size_t theta_steps = 5;

    std::array<vector3, phi_steps * theta_steps> momenta{};

    // Loops of theta values ]0,pi[
    for (std::size_t itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.01 + itheta * (M_PI - 0.01) / theta_steps;

        // Loops of phi values [-pi, pi]
        for (std::size_t iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2. * M_PI) / phi_steps;

            // intialize a track
            vector3 mom{std::cos(phi) * std::sin(theta),
                        std::sin(phi) * std::sin(theta), std::cos(theta)};
            vector::normalize(mom);
            free_track_parameters traj({0., 0., 0.}, 0, mom, -1);

            momenta[itheta * phi_steps + iphi] = traj.mom();
        }
    }

    // Now run the track generator and compare
    std::size_t n_tracks = 0;
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = track.mom();

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Track: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Genrate rays
    n_tracks = 0;
    for (const auto r : uniform_track_generator<ray>(theta_steps, phi_steps)) {
        vector3 &expected = momenta[n_tracks];
        vector3 result = r.dir();

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Ray: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);

    // Generate helix trajectories
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    2. * unit_constants::T};
    n_tracks = 0;
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps)) {
        helix helix_traj(track, &B);
        vector3 &expected = momenta[n_tracks];
        vector3 result = helix_traj.dir(0.);

        EXPECT_NEAR(getter::norm(expected - result), 0, epsilon)
            << "Helix: \n"
            << expected[0] << "\t" << result[0] << "\n"
            << expected[1] << "\t" << result[1] << "\n"
            << expected[2] << "\t" << result[2] << std::endl;

        ++n_tracks;
    }
    ASSERT_EQ(momenta.size(), n_tracks);
}