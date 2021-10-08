/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>

#include "geometry/unified_index_geometry.hpp"
#include "tools/geometry_graph.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the construction of a detector class
TEST(ALGEBRA_PLUGIN, index_geometry) {
    using namespace detray;
    using namespace __plugin;

    using geometry = unified_index_geometry<>;
    using portal = geometry::portal;
    using surface = geometry::surface;
    using graph = geometry_graph<geometry>;

    geometry geo = geometry();

    // Add two volumes
    auto &v0 = geo.new_volume({0., 10., -5., 5., -M_PI, M_PI});
    auto &v1 = geo.new_volume({0., 5., -10., 10., -M_PI, M_PI});

    /// volume 0
    /// volume portals
    auto pt0 = portal(0, {geometry::e_portal_cylinder3, 0}, 0, 1);
    auto pt1 = portal(1, {geometry::e_portal_ring2, 0}, 0, 2);
    auto pt2 = portal(2, {geometry::e_portal_ring2, 1}, 0, 3);
    /// Surface 0
    auto sf0 = surface(3, {geometry::e_rectangle2, 0}, 0, 4);
    /// Surface 1
    auto sf1 = surface(4, {geometry::e_trapezoid2, 0}, 0, 5);

    /// volume 1
    /// volume portals
    auto pt3 = portal(0, {geometry::e_portal_cylinder3, 0}, 1, 6);
    auto pt4 = portal(0, {geometry::e_portal_cylinder3, 1}, 1, 6);
    auto pt5 = portal(1, {geometry::e_portal_ring2, 0}, 1, 7);
    auto pt6 = portal(1, {geometry::e_portal_ring2, 1}, 1, 7);
    /// Surface 2
    auto sf2 = surface(0, {geometry::e_rectangle2, 0}, 1, 8);
    /// Surface 3
    auto sf3 = surface(1, {geometry::e_trapezoid2, 0}, 1, 9);

    /// Correct links for volume 1
    dindex trf_offset_vol1 = 5;
    dindex mask_offset_cyl = 1;
    dindex mask_offset_rg = 2;
    dindex mask_offset_rect = 1;
    dindex mask_offset_trap = 1;

    // update transform links
    geo.update_transform_link(pt3, trf_offset_vol1);
    geo.update_transform_link(pt4, trf_offset_vol1);
    geo.update_transform_link(pt5, trf_offset_vol1);
    geo.update_transform_link(pt6, trf_offset_vol1++);
    geo.update_transform_link(sf2, ++trf_offset_vol1);
    geo.update_transform_link(sf3, trf_offset_vol1);

    // update mask links
    geo.update_mask_link(pt3, mask_offset_cyl);
    geo.update_mask_link(pt4, mask_offset_cyl);
    geo.update_mask_link(pt5, mask_offset_rg);
    geo.update_mask_link(pt6, mask_offset_rg);
    geo.update_mask_link(sf2, mask_offset_rect);
    geo.update_mask_link(sf3, mask_offset_trap);

    /// fill surfaces and portals (through volume alias names)
    dvector<portal> portals_vol0 = {pt0, pt1, pt2};
    dvector<portal> portals_vol1 = {pt3, pt4, pt5, pt6};
    dvector<surface> surfaces_vol0 = {sf0, sf1};
    dvector<surface> surfaces_vol1 = {sf2, sf3};

    auto &v2 = geo.volume_by_index(0);
    auto &v3 = geo.volume_by_index(1);
    geo.template add_objects<geometry::e_portal>(v2, portals_vol0);
    geo.template add_objects<geometry::e_surface>(v2, surfaces_vol0);
    geo.template add_objects<geometry::e_portal>(v3, portals_vol1);
    geo.template add_objects<geometry::e_surface>(v3, surfaces_vol1);

    auto objects_range = darray<dindex, 2>{0, 3};
    ASSERT_TRUE(v2.template range<geometry::e_portal>() == objects_range);

    // Build the graph
    graph g = graph(geo.volumes(), geo.template objects<geometry::e_portal>());

    // Is everything accessible from the graph?
    EXPECT_EQ(g.n_nodes(), geo.n_volumes());
    EXPECT_EQ(g.n_edges(), geo.template n_objects<geometry::e_portal>());

    std::cout << g.to_string() << std::endl;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
