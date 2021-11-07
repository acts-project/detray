/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <unordered_set>

#include "detray/geometry/unified_index_geometry.hpp"
#include "detray/tools/geometry_graph.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the linking of a geometry by loading it into a graph structure
TEST(ALGEBRA_PLUGIN, geometry_linking) {
    using namespace detray;
    using namespace __plugin;

    using geometry = unified_index_geometry<>;
    using portal = geometry::portal;
    using surface = geometry::surface;

    /// Prints linking information for every node when visited
    struct volume_printout {
        void operator()(const geometry::volume_type &n) const {
            std::cout << "On volume: " << n.index() << std::endl;
        }
    };

    using graph = geometry_graph<geometry, volume_printout>;

    geometry geo = geometry();

    // Add two volumes
    geo.new_volume({0., 10., -5., 5., -M_PI, M_PI});
    geo.new_volume({0., 5., -10., 10., -M_PI, M_PI});

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

    // Set the links between volumes {edge, surface_finder}
    // Link to other volume
    pt0.set_edge({1, 1});
    pt1.set_edge({1, 1});
    pt2.set_edge({1, 1});
    // Link back to volume they belong to
    sf0.set_edge({0, 0});
    sf1.set_edge({0, 0});

    pt3.set_edge({0, 0});
    pt4.set_edge({0, 0});
    pt5.set_edge({0, 0});
    pt6.set_edge({0, 0});

    sf2.set_edge({1, 1});
    sf3.set_edge({1, 1});

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
    geo.add_objects(v2, portals_vol0);
    geo.add_objects(v2, surfaces_vol0);
    geo.add_objects(v3, portals_vol1);
    geo.add_objects(v3, surfaces_vol1);

    auto objects_range = darray<dindex, 2>{0, 5};
    ASSERT_TRUE(v2.range() == objects_range);
    objects_range = darray<dindex, 2>{5, 11};
    ASSERT_TRUE(v3.range() == objects_range);

    // Build the graph
    graph g = graph(geo.volumes(), geo.objects());

    // Is everything accessible from the graph?
    EXPECT_EQ(g.n_nodes(), geo.n_volumes());
    EXPECT_EQ(g.n_edges(), geo.n_objects());

    std::cout << g.to_string() << std::endl;
    std::cout << "Walking through geometry: " << std::endl;
    g.bfs();

    const auto &adj = g.adjacency_list();

    // Volume 0 has 3 portals to volume 1 and two surfaces linking to itself
    std::map<dindex, dindex> nbrs_map_v0{{0, 2}, {1, 3}};
    // Volume 1 has 4 portals to volume 0 and two surfaces linking to itself
    std::map<dindex, dindex> nbrs_map_v1{{0, 4}, {1, 2}};

    // Check this with graph
    ASSERT_TRUE(adj.at(0) == nbrs_map_v0);
    ASSERT_TRUE(adj.at(1) == nbrs_map_v1);
    // ASSERT_TRUE(adj.at(0).first == obj_hash_v0);
    // ASSERT_TRUE(adj.at(1).first == obj_hash_v1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
