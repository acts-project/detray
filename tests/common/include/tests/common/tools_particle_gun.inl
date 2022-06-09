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
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/test_trajectories.hpp"
#include "tests/common/tools/track_generators.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

constexpr const float epsilon = 1e-3;

/// Re-use the intersection kernel test for particle gun
TEST(tools, helix_intersector) {

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
    using mask_container_t = typename mask_defs::template mask_store_type<>;

    /// The Surface definition:
    /// <transform_link, volume_link, source_link, link_type_in_mask>
    using surface_t = surface<mask_defs, dindex, dindex, source_link_t>;
    using surface_container_t = dvector<surface_t>;

    // The transforms & their store
    static_transform_store<>::context static_context{};
    static_transform_store transform_store;
    // Transforms of the rectangle, trapezoid and annulus
    transform_store.emplace_back(static_context, point3{0., 0., 10.});
    transform_store.emplace_back(static_context, point3{0., 0., 20.});
    transform_store.emplace_back(static_context, point3{0., -20., 30.});
    // The masks & their store
    mask_container_t mask_store(host_mr);
    mask_store.template add_value<e_rectangle2>(10., 10., 0);
    mask_store.template add_value<e_trapezoid2>(10., 20., 30., 0);
    mask_store.template add_value<e_annulus2>(15., 55., 0.75, 1.95, 2., -2., 0.,
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
    helix h({pos, 0, mom, -1}, &B);

    // Validation data
    point3 expected_rectangle{0.01, 0.01, 10.};
    point3 expected_trapezoid{0.02, 0.02, 20.};
    point3 expected_annulus{0.03, 0.03, 30.};
    std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersections - with automated dispatching via the kernel
    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        auto sfi_helix = particle_gun::intersect(h, surface, transform_store,
                                                 mask_store, epsilon);

        ASSERT_NEAR(sfi_helix.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[2], expected_points[sf_idx][2], 1e-7);

        auto sfi_ray = particle_gun::intersect(r, surface, transform_store,
                                               mask_store, epsilon);
        ASSERT_NEAR(sfi_ray.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi_ray.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi_ray.p3[2], expected_points[sf_idx][2], 1e-7);
    }
}

/// Brute force test: Intersect toy geometry and compare between ray and helix
/// without B-field
TEST(tools, particle_gun) {

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);

    unsigned int theta_steps = 50;
    unsigned int phi_steps = 50;
    const point3 ori{0., 0., 0.};

    // Record ray tracing
    std::vector<std::vector<std::pair<dindex, particle_gun::intersection_type>>>
        expected;
    //  Iterate through uniformly distributed momentum directions with ray
    for (const auto test_ray :
         uniform_track_generator<ray>(theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            particle_gun::shoot_particle(toy_det, test_ray, epsilon);

        expected.push_back(intersection_record);
    }

    // Simulate straight line track
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    0.00001 * unit_constants::T};
    // Iterate through uniformly distributed momentum directions with helix
    std::size_t n_tracks{0};
    for (const auto track : uniform_track_generator<free_track_parameters>(
             theta_steps, phi_steps, ori)) {
        helix test_helix(track, &B);

        // Record all intersections and objects along the ray
        const auto intersection_trace =
            particle_gun::shoot_particle(toy_det, test_helix, epsilon);

        // Should have encountered the same number of tracks (vulnerable to
        // floating point errors)
        EXPECT_EQ(expected[n_tracks].size(), intersection_trace.size());

        // Check every single recorded intersection
        for (std::size_t i = 0; i < intersection_trace.size(); ++i) {
            if (expected[n_tracks][i].first != intersection_trace[i].first) {
                // Intersection record at portal bound might be flipped
                // (the portals overlap completely)
                if (expected[n_tracks][i].first ==
                        intersection_trace[i + 1].first and
                    expected[n_tracks][i + 1].first ==
                        intersection_trace[i].first) {
                    // Have already checked the next record
                    ++i;
                    continue;
                }
            }
            EXPECT_EQ(expected[n_tracks][i].first, intersection_trace[i].first);
        }

        ++n_tracks;
    }
}