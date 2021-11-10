/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/core/mask_store.hpp"
#include "detray/core/track.hpp"
#include "detray/core/transform_store.hpp"
#include "detray/geometry/surface_base.hpp"
#include "detray/masks/masks.hpp"
#include "detray/tools/concentric_cylinder_intersector.hpp"
#include "detray/tools/cylinder_intersector.hpp"
#include "detray/tools/intersection_kernel.hpp"
#include "detray/tools/planar_intersector.hpp"
#include "detray/utils/enumerate.hpp"

using namespace detray;
using namespace __plugin;

// This tests the construction of a surface
TEST(tools, intersection_kernel_single) {

    vecmem::host_memory_resource host_mr;

    /// Surface components:
    using mask_link = darray<dindex, 1>;
    using surface_link = dindex;
    /// - masks, with mask identifiers 0,1,2
    using surface_rectangle =
        rectangle2<planar_intersector, __plugin::cartesian2, mask_link, 0>;
    using surface_trapezoid =
        trapezoid2<planar_intersector, __plugin::cartesian2, mask_link, 1>;
    using surface_annulus =
        annulus2<planar_intersector, __plugin::cartesian2, mask_link, 2>;

    /// - mask index: type, entry
    using surface_mask_index = darray<dindex, 2>;
    using surface_mask_container =
        mask_store<dtuple, dvector, surface_rectangle, surface_trapezoid,
                   surface_annulus>;

    /// The Surface definition:
    /// <transform_link, mask_link, volume_link, source_link, link_type_in_mask>
    using surface = surface_base<dindex, surface_mask_index, dindex,
                                 surface_link, mask_link>;
    using surface_container = dvector<surface>;

    // The transforms & their store
    transform3 rectangle_transform(point3{0., 0., 10.});
    transform3 trapezoid_transform(point3{0., 0., 20.});
    transform3 annulus_transform(point3{0., -20., 30.});
    static_transform_store<>::context static_context;
    static_transform_store transform_store;
    transform_store.push_back(static_context, rectangle_transform);
    transform_store.push_back(static_context, trapezoid_transform);
    transform_store.push_back(static_context, annulus_transform);
    // The masks & their store
    surface_mask_container mask_store(host_mr);
    mask_store.template add_mask<0>(10., 10.);
    mask_store.template add_mask<1>(10., 20., 30.);
    mask_store.template add_mask<2>(15., 55., 0.75, 1.95, 2., -2.);
    // The surfaces and their store
    surface rectangle_surface(0u, {0, 0}, 0, 0);
    surface trapezoid_surface(1u, {1, 0}, 0, 1);
    surface annulus_surface(2u, {2, 0}, 0, 2);
    surface_container surfaces = {rectangle_surface, trapezoid_surface,
                                  annulus_surface};

    track<decltype(transform_store)::context> track;
    track.pos = point3{0., 0., 0.};
    track.dir = vector::normalize(vector3{0.01, 0.01, 10.});

    // Validation data
    point3 expected_rectangle{0.01, 0.01, 10.};
    point3 expected_trapezoid{0.02, 0.02, 20.};
    point3 expected_annulus{0.03, 0.03, 30.};

    std::vector<point3> expected_points = {
        expected_rectangle, expected_trapezoid, expected_annulus};

    // Try the intersection - with automated dispatching via the kernel
    unsigned int it = 0;
    for (const auto &surface : surfaces) {
        auto sfi = intersect(track, surface, transform_store, mask_store);

        ASSERT_NEAR(sfi.p3[0], expected_points[it][0], 1e-7);
        ASSERT_NEAR(sfi.p3[1], expected_points[it][1], 1e-7);
        ASSERT_NEAR(sfi.p3[2], expected_points[it][2], 1e-7);
        ++it;
    }

    return;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
