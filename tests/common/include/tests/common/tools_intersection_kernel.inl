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
#include "detray/intersection/concentric_cylinder_intersector.hpp"
#include "detray/intersection/cylinder_intersector.hpp"
#include "detray/intersection/intersection_kernel.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/propagator/track.hpp"
#include "detray/utils/enumerate.hpp"

using namespace detray;

// This tests the construction of a surface
TEST(tools, intersection_kernel_single) {

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
    free_track_parameters track(pos, 0, mom, -1);

    // Validation data
    point3 expected_rectangle{0.01, 0.01, 10.};
    point3 expected_trapezoid{0.02, 0.02, 20.};
    point3 expected_annulus{0.03, 0.03, 30.};

    std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersection - with automated dispatching via the kernel
    unsigned int it = 0;
    for (const auto &_surface : surfaces) {
        auto sfi =
            intersection_kernel{}(_surface, track, transform_store, mask_store);

        ASSERT_NEAR(sfi.p3[0], expected_points[it][0], 1e-7);
        ASSERT_NEAR(sfi.p3[1], expected_points[it][1], 1e-7);
        ASSERT_NEAR(sfi.p3[2], expected_points[it][2], 1e-7);
        ++it;
    }

    return;
}
