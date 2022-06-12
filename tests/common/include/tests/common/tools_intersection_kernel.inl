/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/mask_store.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/core/type_registry.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/intersection/helix_plane_intersector.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/ray_plane_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/enumerate.hpp"

using namespace detray;

constexpr const float epsilon = 1e-3;

// This tests the construction of a surface
TEST(tools, intersection_kernel_ray) {

    using vector3 = __plugin::vector3<scalar>;
    using point3 = __plugin::point3<scalar>;

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
        rectangle2<ray_plane_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_t>;
    using trapezoid_t =
        trapezoid2<ray_plane_intersector, __plugin::cartesian2<detray::scalar>,
                   edge_t>;
    using annulus_t = annulus2<ray_plane_intersector,
                               __plugin::cartesian2<detray::scalar>, edge_t>;

    using mask_defs =
        mask_registry<mask_ids, rectangle_t, trapezoid_t, annulus_t>;
    using mask_container_t = typename mask_defs::mask_store_type<>;

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
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, 0, 0, false);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, 0, 1, false);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, 0, 2, false);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};

    const point3 pos{0., 0., 0.};
    const vector3 mom{0.01, 0.01, 10.};
    const free_track_parameters track(pos, 0, mom, -1);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};

    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersection - with automated dispatching via the kernel
    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        const auto sfi = intersect(track, surface, transform_store, mask_store);

        ASSERT_NEAR(sfi.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi.p3[2], expected_points[sf_idx][2], 1e-7);
    }
}

/// Re-use the intersection kernel test for particle gun
TEST(tools, intersection_kernel_helix) {

    using vector3 = __plugin::vector3<scalar>;
    using point3 = __plugin::point3<scalar>;

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
        rectangle2<helix_plane_intersector,
                   __plugin::cartesian2<detray::scalar>, edge_t>;
    using trapezoid_t =
        trapezoid2<helix_plane_intersector,
                   __plugin::cartesian2<detray::scalar>, edge_t>;
    using annulus_t = annulus2<helix_plane_intersector,
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
    const surface_t rectangle_surface(0u, {e_rectangle2, 0}, 0, 0, false);
    const surface_t trapezoid_surface(1u, {e_trapezoid2, 0}, 0, 1, false);
    const surface_t annulus_surface(2u, {e_annulus2, 0}, 0, 2, false);
    surface_container_t surfaces = {rectangle_surface, trapezoid_surface,
                                    annulus_surface};
    const point3 pos{0., 0., 0.};
    const vector3 mom{0.01, 0.01, 10.};
    const vector3 B{0. * unit_constants::T, 0. * unit_constants::T,
                    epsilon * unit_constants::T};
    const detail::helix h({pos, 0, mom, -1}, &B);

    // Validation data
    const point3 expected_rectangle{0.01, 0.01, 10.};
    const point3 expected_trapezoid{0.02, 0.02, 20.};
    const point3 expected_annulus{0.03, 0.03, 30.};
    const std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersections - with automated dispatching via the kernel
    for (const auto& [sf_idx, surface] : enumerate(surfaces)) {
        const auto sfi_helix =
            intersect(h, surface, transform_store, mask_store);

        ASSERT_NEAR(sfi_helix.p3[0], expected_points[sf_idx][0], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[1], expected_points[sf_idx][1], 1e-7);
        ASSERT_NEAR(sfi_helix.p3[2], expected_points[sf_idx][2], 1e-7);
    }
}
