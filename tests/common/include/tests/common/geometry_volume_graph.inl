/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/tools/volume_graph.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the linking of a geometry by loading it into a graph structure
TEST(ALGEBRA_PLUGIN, geometry_linking) {
    using namespace detray;
    using namespace __plugin;

    vecmem::host_memory_resource host_mr;
    auto det = create_toy_geometry(host_mr);
    using detector_t = decltype(det);

    /// Prints linking information for every node when visited
    struct volume_printout {
        void operator()(const detector_t::volume_type &n) const {
            std::cout << "On volume: " << n.index() << std::endl;
        }
    };

    using graph = volume_graph<detector_t, volume_printout>;

    // Build the graph
    graph g = graph(det.volumes(), det.surfaces(), det.mask_store());

    // Is everything accessible from the graph?
    EXPECT_EQ(g.n_nodes(), det.volumes().size());

    std::cout << g.to_string() << std::endl;
    std::cout << "Walking through geometry: " << std::endl;
    // g.bfs();

    // Volume 0 has 3 portals to volume 1 and two surfaces linking to itself
    // std::map<dindex, dindex> nbrs_map_v0{{0, 2}, {1, 3}};
    // Volume 1 has 4 portals to volume 0 and two surfaces linking to itself
    // std::map<dindex, dindex> nbrs_map_v1{{0, 4}, {1, 2}};

    // Check this with graph
    // const auto &adj = g.adjacency_list();
    // ASSERT_TRUE(adj.at(0) == nbrs_map_v0);
    // ASSERT_TRUE(adj.at(1) == nbrs_map_v1);
}
