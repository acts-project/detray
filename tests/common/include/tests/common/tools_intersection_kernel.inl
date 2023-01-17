/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/cylinder_portal_intersector.hpp"
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

enum mask_ids : unsigned int {
    e_rectangle2 = 0,
    e_trapezoid2 = 1,
    e_annulus2 = 2,
    e_cylinder2 = 3,
    e_cylinder2_portal = 4,
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
using cylinder_t = mask<cylinder2D<true, cylinder_intersector>, volume_link_t>;
using cylinder_portal_t =
    mask<cylinder2D<false, cylinder_portal_intersector>, volume_link_t>;

using mask_container_t =
    regular_multi_store<mask_ids, empty_context, dtuple, dvector, rectangle_t,
                        trapezoid_t, annulus_t, cylinder_t, cylinder_portal_t>;
using mask_link_t = typename mask_container_t::single_link;
using material_link_t = dtyped_index<material_ids, dindex>;

using transform_container_t = single_store<__plugin::transform3<scalar>>;
using transform3_t = typename transform_container_t::value_type;
using transform_link_t = dindex;

/// The Surface definition:
using surface_t = surface<mask_link_t, material_link_t, transform_link_t>;
using surface_container_t = dvector<surface_t>;

using vector3 = typename transform3_t::vector3;
using point3 = typename transform3_t::point3;

constexpr const scalar epsilon{1e-3f};
constexpr const scalar is_close{1e-6f};

// TODO: How about merging ray and helix tests into one to remove the code
// repetition?

// This tests the construction of a surface
TEST(tools, intersection_kernel_ray) {
    vecmem::host_memory_resource host_mr;

    // The transforms & their store
    typename transform_container_t::context_type static_context{};
    transform_container_t transform_store;
    // Transforms of the rectangle, trapezoid, annulus and the two cylinders
    transform_store.emplace_back(static_context, vector3{0.f, 0.f, 10.f});
    transform_store.emplace_back(static_context, vector3{0.f, 0.f, 20.f});
    transform_store.emplace_back(static_context, vector3{0.f, -20.f, 30.f});
    // 90deg rotation around y-axis
    transform_store.emplace_back(static_context, vector3{0.f, 0.f, 50.f},
                                 vector3{1.f, 0.f, 0.f},
                                 vector3{0.f, 0.f, -1.f});
    transform_store.emplace_back(static_context, vector3{0.f, 0.f, 100.f},
                                 vector3{1.f, 0.f, 0.f},
                                 vector3{0.f, 0.f, -1.f});

    // The masks & their store
    mask_container_t mask_store(host_mr);
    mask_store.template emplace_back<e_rectangle2>(empty_context{}, 0UL, 10.f,
                                                   10.f);
    mask_store.template emplace_back<e_trapezoid2>(empty_context{}, 0UL, 10.f,
                                                   20.f, 30.f);
    mask_store.template emplace_back<e_annulus2>(
        empty_context{}, 0UL, 15.f, 55.f, 0.75f, 1.95f, 2.f, -2.f, 0.f);
    mask_store.template emplace_back<e_cylinder2>(empty_context{}, 0UL, 5.f,
                                                  -10.f, 10.f);
    mask_store.template emplace_back<e_cylinder2_portal>(empty_context{}, 0UL,
                                                         4.f, -10.f, 10.f);

    // The surfaces and their store
    surface_t rectangle_surface(0u, {e_rectangle2, 0}, {e_slab, 0}, 0, 0,
                                surface_id::e_sensitive);
    rectangle_surface.set_barcode(0UL);
    surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, {e_slab, 1}, 0, 1,
                                surface_id::e_sensitive);
    trapezoid_surface.set_barcode(1UL);
    surface_t annulus_surface(2u, {e_annulus2, 0}, {e_slab, 2}, 0, 2,
                              surface_id::e_sensitive);
    annulus_surface.set_barcode(2UL);
    surface_t cyl_surface(3u, {e_cylinder2, 0}, {e_slab, 2}, 0, 3,
                          surface_id::e_passive);
    cyl_surface.set_barcode(3UL);
    surface_t cyl_portal_surface(4u, {e_cylinder2_portal, 0}, {e_slab, 2}, 0, 4,
                                 surface_id::e_portal);
    cyl_portal_surface.set_barcode(4UL);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface, cyl_surface,
                                    cyl_portal_surface};

    const point3 pos{0.f, 0.f, 0.f};
    const vector3 mom{0.01f, 0.01f, 10.f};
    const free_track_parameters<transform3_t> track(pos, 0, mom, -1);

    // Validation data (no valid intersection on the cylinders)
    const point3 expected_rectangle{0.01f, 0.01f, 10.f};
    const point3 expected_trapezoid{0.02f, 0.02f, 20.f};
    const point3 expected_annulus{0.03f, 0.03f, 30.f};
    const point3 expected_cylinder1{0.045f, 0.045f, 45.000195f};
    const point3 expected_cylinder2{0.055f, 0.055f, 54.999706f};
    const point3 expected_cylinder_pt{0.096001f, 0.096001f, 96.0011215f};

    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus,
        expected_cylinder1, expected_cylinder2, expected_cylinder_pt};

    // Initialize kernel
    std::vector<line_plane_intersection<surface_t, transform3_t>> sfi_init;

    for (const auto &surface : surfaces) {
        mask_store.visit<intersection_initialize>(surface.mask(), sfi_init,
                                                  detail::ray(track), surface,
                                                  transform_store, epsilon);
    }
    // Also check intersections
    for (std::size_t i = 0UL; i < expected_points.size(); ++i) {
        ASSERT_EQ(sfi_init[i].direction, intersection::direction::e_along);
        ASSERT_EQ(sfi_init[i].volume_link, 0UL);
        ASSERT_NEAR(sfi_init[i].p3[0], expected_points[i][0], is_close)
            << " at surface " << sfi_init[i].surface.barcode();
        ASSERT_NEAR(sfi_init[i].p3[1], expected_points[i][1], is_close)
            << " at surface " << sfi_init[i].surface.barcode();
        ASSERT_NEAR(sfi_init[i].p3[2], expected_points[i][2], is_close)
            << " at surface " << sfi_init[i].surface.barcode();
    }

    // Update kernel
    /*std::vector<line_plane_intersection<surface_t, transform3_t>> sfi_update;
    sfi_update.resize(5);

    for (const auto [idx, surface] : detray::views::enumerate(surfaces)) {
        sfi_update[idx].surface = surface;
        mask_store.visit<intersection_update>(
            surface.mask(), detail::ray(track), sfi_update[idx],
    transform_store);

        if(sfi_init[idx].status != intersection::status::e_inside) { continue; }
        ASSERT_EQ(sfi_update[idx].direction, intersection::direction::e_along);
        ASSERT_EQ(sfi_update[idx].volume_link, 0UL);
        ASSERT_NEAR(sfi_update[idx].p3[0], expected_points[idx][0], is_close) <<
    " at surface " << idx;; ASSERT_NEAR(sfi_update[idx].p3[1],
    expected_points[idx][1], is_close) << " at surface " << idx;;
        ASSERT_NEAR(sfi_update[idx].p3[2], expected_points[idx][2], is_close) <<
    " at surface " << idx;;
    }

    // Compare
    ASSERT_EQ(sfi_init.size(), 5);
    ASSERT_EQ(sfi_update.size(), 5);
    for (int i = 0; i < 5; i++) {
        ASSERT_EQ(sfi_init[i].p3, sfi_update[i].p3);
        ASSERT_EQ(sfi_init[i].p2, sfi_update[i].p2);
        ASSERT_EQ(sfi_init[i].path, sfi_update[i].path);
    }*/
}

/// Re-use the intersection kernel test for particle gun
TEST(tools, intersection_kernel_helix) {

    vecmem::host_memory_resource host_mr;

    // The transforms & their store
    typename transform_container_t::context_type static_context{};
    transform_container_t transform_store;
    // Transforms of the rectangle, trapezoid and annulus
    transform_store.emplace_back(static_context, vector3{0., 0., 10.});
    transform_store.emplace_back(static_context, vector3{0., 0., 20.});
    transform_store.emplace_back(static_context, vector3{0., -20., 30.});
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
    const detail::helix<transform3_t> h({pos, 0, mom, -1}, &B);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};
    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};
    std::vector<line_plane_intersection<surface_t, transform3_t>> sfi_helix{};

    // Try the intersections - with automated dispatching via the kernel
    for (const auto [sf_idx, surface] : detray::views::enumerate(surfaces)) {
        mask_store.visit<helix_intersection_initialize>(
            surface.mask(), sfi_helix, h, surface, transform_store);

        ASSERT_NEAR(sfi_helix[0].p3[0], expected_points[sf_idx][0], is_close);
        ASSERT_NEAR(sfi_helix[0].p3[1], expected_points[sf_idx][1], is_close);
        ASSERT_NEAR(sfi_helix[0].p3[2], expected_points[sf_idx][2], is_close);

        sfi_helix.clear();
    }
}
