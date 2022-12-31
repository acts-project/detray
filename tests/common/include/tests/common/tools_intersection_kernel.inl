/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/ranges.hpp"
#include "tests/common/tools/intersectors/helix_intersection_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Test include(s)
#include <gtest/gtest.h>

using namespace detray;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3_type = __plugin::transform3<scalar>;

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_trapezoid2 = 1,
    e_annulus2 = 2,
};

enum material_ids : unsigned int {
    e_slab = 0,
};

/// Surface components:
using volume_link_t = dindex;

/// - masks, with mask identifiers 0,1,2
using rectangle_t = mask<rectangle2D<>, volume_link_t>;
using trapezoid_t = mask<trapezoid2D<>, volume_link_t>;
using annulus_t = mask<annulus2D<>, volume_link_t>;

using mask_container_t =
    regular_multi_store<mask_ids, empty_context, dtuple, dvector, rectangle_t,
                        trapezoid_t, annulus_t>;
using mask_link_t = typename mask_container_t::single_link;
using material_link_t = dtyped_index<material_ids, dindex>;

using transform_container_t = single_store<__plugin::transform3<scalar>>;
using transform_link_t = dindex;

/// The Surface definition:
using surface_t = surface<mask_link_t, material_link_t, transform_link_t>;
using surface_container_t = dvector<surface_t>;

constexpr const float epsilon = 1e-3;

// TODO: How about merging ray and helix tests into one to remove the code
// repetition?

// This tests the construction of a surface
TEST(tools, intersection_kernel_ray) {
    vecmem::host_memory_resource host_mr;

    // The transforms & their store
    typename transform_container_t::context_type static_context{};
    transform_container_t transform_store;
    // Transforms of the rectangle, trapezoid and annulus
    transform_store.emplace_back(static_context, point3{0., 0., 10.});
    transform_store.emplace_back(static_context, point3{0., 0., 20.});
    transform_store.emplace_back(static_context, point3{0., -20., 30.});
    // The masks & their store
    mask_container_t mask_store(host_mr);
    mask_store.template emplace_back<e_rectangle2>(empty_context{}, 0UL, 10.f,
                                                   10.f);
    mask_store.template emplace_back<e_trapezoid2>(empty_context{}, 0UL, 10.f,
                                                   20.f, 30.f);
    mask_store.template emplace_back<e_annulus2>(
        empty_context{}, 0UL, 15.f, 55.f, 0.75f, 1.95f, 2.f, -2.f, 0.f);

    // The surfaces and their store
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, {e_slab, 0}, 0, 0,
                                      surface_id::e_sensitive);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, {e_slab, 1}, 0, 1,
                                      surface_id::e_sensitive);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, {e_slab, 2}, 0, 2,
                                    surface_id::e_sensitive);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};

    const point3 pos{0., 0., 0.};
    const vector3 mom{0.01, 0.01, 10.};
    const free_track_parameters<transform3_type> track(pos, 0, mom, -1);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};

    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Initialize kernel
    std::vector<line_plane_intersection> sfi_init;

    for (const auto& surface : surfaces) {
        mask_store.visit<intersection_initialize>(surface.mask(), sfi_init,
                                                  detail::ray(track), surface,
                                                  transform_store);
    }

    // Update kernel
    std::vector<line_plane_intersection> sfi_update;

    for (const auto [sf_idx, surface] : detray::views::enumerate(surfaces)) {
        const auto sfi = mask_store.visit<intersection_update>(
            surface.mask(), detail::ray(track), surface, transform_store);

        sfi_update.push_back(sfi);

        ASSERT_NEAR(sfi.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi.p3[2], expected_points[sf_idx][2], 1e-7);
    }

    // Compare
    ASSERT_EQ(sfi_init.size(), 3);
    ASSERT_EQ(sfi_update.size(), 3);
    for (int i = 0; i < 3; i++) {
        ASSERT_EQ(sfi_init[i].p3, sfi_update[i].p3);
        ASSERT_EQ(sfi_init[i].p2, sfi_update[i].p2);
        ASSERT_EQ(sfi_init[i].path, sfi_update[i].path);
    }
}

/// Re-use the intersection kernel test for particle gun
TEST(tools, intersection_kernel_helix) {

    vecmem::host_memory_resource host_mr;

    // The transforms & their store
    typename transform_container_t::context_type static_context{};
    transform_container_t transform_store;
    // Transforms of the rectangle, trapezoid and annulus
    transform_store.emplace_back(static_context, point3{0., 0., 10.});
    transform_store.emplace_back(static_context, point3{0., 0., 20.});
    transform_store.emplace_back(static_context, point3{0., -20., 30.});
    // The masks & their store
    mask_container_t mask_store(host_mr);
    mask_store.template emplace_back<e_rectangle2>(empty_context{}, 0UL, 10.f,
                                                   10.f);
    mask_store.template emplace_back<e_trapezoid2>(empty_context{}, 0UL, 10.f,
                                                   20.f, 30.f);
    mask_store.template emplace_back<e_annulus2>(
        empty_context{}, 0UL, 15.f, 55.f, 0.75f, 1.95f, 2.f, -2.f, 0.f);

    // The surfaces and their store
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, {e_slab, 0}, 0, 0,
                                      surface_id::e_sensitive);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, {e_slab, 1}, 0, 1,
                                      surface_id::e_sensitive);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, {e_slab, 2}, 0, 2,
                                    surface_id::e_sensitive);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};
    const point3 pos{0., 0., 0.};
    const vector3 mom{0.01, 0.01, 10.};
    const vector3 B{0. * unit<scalar>::T, 0. * unit<scalar>::T,
                    epsilon * unit<scalar>::T};
    const detail::helix<transform3_type> h({pos, 0, mom, -1}, &B);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};
    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersections - with automated dispatching via the kernel
    for (const auto [sf_idx, surface] : detray::views::enumerate(surfaces)) {
        const auto sfi_helix = mask_store.visit<helix_intersection_update>(
            surface.mask(), h, surface, transform_store);

        ASSERT_NEAR(sfi_helix.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[2], expected_points[sf_idx][2], 1e-7);
    }
}
