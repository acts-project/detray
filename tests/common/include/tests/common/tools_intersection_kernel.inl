/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2021 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "core/surface_base.hpp"
#include "core/track.hpp"
#include "core/transform_store.hpp"
#include "masks/masks.hpp"
#include "tools/intersection_kernel.hpp"
#include "tools/planar_intersector.hpp"
#include "tools/cylinder_intersector.hpp"
#include "tools/concentric_cylinder_intersector.hpp"
#include "utils/enumerate.hpp"

#include <gtest/gtest.h>

using namespace detray;
using namespace __plugin;

// This tests the construction of a surface
TEST(tools, intersection_kernel_single)
{

    using vector3 = transform3::vector3;
    using point3 = transform3::point3;
    using point2 = cartesian2::point2;

    /// Surface components:
    using mask_link = dindex;
    using surface_link = dindex;
    /// - masks, with mask identifiers 0,1,2
    using surface_rectangle = rectangle2<planar_intersector, __plugin::cartesian2, mask_link, 0>;
    using surface_trapezoid = trapezoid2<planar_intersector, __plugin::cartesian2, mask_link, 1>;
    using surface_annulus = annulus2<planar_intersector, __plugin::cartesian2, mask_link, 2>;
    using surface_cylinder = cylinder3<false, cylinder_intersector, __plugin::cylindrical2, mask_link, 3>;
    using surface_concentric_cylinder = cylinder3<false, concentric_cylinder_intersector, __plugin::cylindrical2, mask_link, 4>;
    /// - mask index: type, entry
    using surface_mask_index = darray<dindex, 2>;
    using surface_mask_container = dtuple<dvector<surface_rectangle>, 
                                 dvector<surface_trapezoid>, 
                                 dvector<surface_annulus>, 
                                 dvector<surface_cylinder>,
                                 dvector<surface_concentric_cylinder> >;

    /// The Surface definition:
    ///  <transform_link, mask_link, volume_link, source_link >
    using surface = surface_base<dindex, surface_mask_index, dindex, surface_link>;
    using surface_container = dvector<surface>;

    using surface_intersection = intersection<point3, point2>;

    // The transforms & their store
    transform3 rectangle_transform(point3{0., 0., 10.});
    transform3 trapezoid_transform(point3{0., 0., 20.});
    transform3 annulus_transform(point3{0., -20., 30.});
    static_transform_store::context static_context;
    static_transform_store transform_store;
    transform_store.push_back(static_context, rectangle_transform);
    transform_store.push_back(static_context, trapezoid_transform);
    transform_store.push_back(static_context, annulus_transform);
    // The masks & their store
    surface_rectangle first_rectangle = {10., 10.};
    surface_trapezoid second_trapezoid = {10., 20., 30.};
    surface_annulus thrid_annulus = {15., 55., 0.75, 1.95, 2., -2.};
    surface_mask_container mask_store;
    std::get<0>(mask_store).push_back(first_rectangle);
    std::get<1>(mask_store).push_back(second_trapezoid);
    std::get<2>(mask_store).push_back(thrid_annulus);
    // The surfaces and their store
    surface rectangle_surface(0u, {0, 0}, 0, 0);
    surface trapezoid_surface(1u, {1, 0}, 0, 1);
    surface annulus_surface(2u, {2, 0}, 0, 2);
    surface_container surfaces = {rectangle_surface, trapezoid_surface, annulus_surface};

    // Try the intersection - first one by one

    track<transform3, decltype(transform_store)::context> track;
    track.pos = point3{0., 0., 0.};
    track.dir = vector::normalize(vector3{0.01, 0.01, 10.});

    // Quick helper to check for within epsilon
    auto within_epsilon = [](const point3 &a, const point3 &b, scalar epsilon) -> bool {
        return (std::abs(a[0] - b[0]) < epsilon && std::abs(a[1] - b[1]) < epsilon && std::abs(a[2] - b[2]) < epsilon);
    };

    // Intersect the first surface
    auto sfi_rectangle = intersect_by_group(
        track, rectangle_transform, std::get<0>(mask_store), 0);

    point3 expected_rectangle{0.01, 0.01, 10.};
    ASSERT_TRUE(within_epsilon(std::get<0>(sfi_rectangle).point3, expected_rectangle, 1e-7));

    auto sfi_trapezoid = intersect_by_group(
        track, trapezoid_transform, std::get<1>(mask_store), 0);

    point3 expected_trapezoid{0.02, 0.02, 20.};
    ASSERT_TRUE(within_epsilon(std::get<0>(sfi_trapezoid).point3, expected_trapezoid, 1e-7));

    auto sfi_annulus = intersect_by_group(
        track, annulus_transform, std::get<2>(mask_store), 0);

    point3 expected_annulus{0.03, 0.03, 30.};
    ASSERT_TRUE(within_epsilon(std::get<0>(sfi_annulus).point3, expected_annulus, 1e-7));

    std::vector<point3> expected_points = {expected_rectangle, expected_trapezoid, expected_annulus};
    std::vector<point3> result_points = {};

    // Try the intersection - with automated dispatching via the kernel
    unsigned int it = 0;
    for (const auto &surface : surfaces)
    {
        auto sfi_surface = intersect(track, surface, transform_store, mask_store);

        const auto &sfi = std::get<0>(sfi_surface);
        result_points.push_back(sfi.point3);

        ASSERT_TRUE(within_epsilon(result_points[it], expected_points[it], 1e-7));
        ++it;
    }

    return;
}


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
