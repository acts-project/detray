/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/enumerate.hpp"
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
using edge_t = dindex;
using source_link_t = dindex;

/// - masks, with mask identifiers 0,1,2
using rectangle_t =
    rectangle2<transform3_type, plane_intersector, cartesian2, edge_t>;
using trapezoid_t =
    trapezoid2<transform3_type, plane_intersector, cartesian2, edge_t>;
using annulus_t =
    annulus2<transform3_type, plane_intersector, cartesian2, edge_t>;

using mask_defs =
    tuple_vector_registry<mask_ids, rectangle_t, trapezoid_t, annulus_t>;
using mask_container_t = typename mask_defs::template store_type<>;

constexpr const float epsilon = 1e-3;

// TODO: How about merging ray and helix tests into one to remove the code
// repetition?

// This tests the construction of a surface
TEST(tools, intersection_kernel_ray) {
    vecmem::host_memory_resource host_mr;

    // Materials with a slab
    using material_defs =
        tuple_vector_registry<material_ids, material_slab<scalar>>;
    using material_container_t = typename material_defs::store_type<>;

    /// The Surface definition:
    /// <transform_link, volume_link, source_link, link_type_in_mask>
    using surface_t = surface<mask_defs, material_defs, dindex, source_link_t>;
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
    // Materials and their store
    material_container_t material_store(host_mr);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              1. * unit_constants::mm);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              2. * unit_constants::mm);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              3. * unit_constants::mm);

    // The surfaces and their store
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, {e_slab, 0}, 0, 0,
                                      false);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, {e_slab, 1}, 0, 1,
                                      false);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, {e_slab, 2}, 0, 2,
                                    false);
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

    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        mask_store.call<intersection_initialize>(surface.mask(), sfi_init,
                                                 detail::ray(track), surface,
                                                 transform_store);
    }

    // Update kernel
    std::vector<line_plane_intersection> sfi_update;

    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        const auto sfi = mask_store.call<intersection_update>(
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

    // Materials with a slab
    using material_defs =
        tuple_vector_registry<material_ids, material_slab<scalar>>;
    using material_container_t = typename material_defs::store_type<>;

    /// The Surface definition:
    /// <transform_link, volume_link, source_link, link_type_in_mask>
    using surface_t = surface<mask_defs, material_defs, dindex, source_link_t>;
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
    // Materials and their store
    material_container_t material_store(host_mr);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              1. * unit_constants::mm);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              2. * unit_constants::mm);
    material_store.template add_value<e_slab>(silicon<scalar>(),
                                              3. * unit_constants::mm);

    // The surfaces and their store
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, {e_slab, 0}, 0, 0,
                                      false);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, {e_slab, 1}, 0, 1,
                                      false);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, {e_slab, 2}, 0, 2,
                                    false);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};
    const point3 pos{0., 0., 0.};
    const vector3 mom{0.01, 0.01, 10.};
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    epsilon * unit_constants::T};
    const detail::helix<transform3_type> h({pos, 0, mom, -1}, &B);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};
    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersections - with automated dispatching via the kernel
    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        const auto sfi_helix = mask_store.call<helix_intersection_update>(
            surface.mask(), h, surface, transform_store);

        ASSERT_NEAR(sfi_helix.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[2], expected_points[sf_idx][2], 1e-7);
    }
}
