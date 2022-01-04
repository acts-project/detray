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

    // Build the graph
    volume_graph<detector_t, volume_printout> graph(det);

    // Is everything accessible from the graph?
    EXPECT_EQ(graph.n_nodes(), det.volumes().size());

    std::cout << graph.to_string() << std::endl;
    std::cout << "Walking through geometry: " << std::endl;
    // graph.bfs();

    const auto &adj_mat = graph.adjacency_matrix();

    // Volume 0 has 3 portals to volume 1 and two surfaces linking to itself
    // Volume 1 has 4 portals to volume 0 and two surfaces linking to itself
    dvector<dindex> adj_mat_truth = {2, 3, 0, 4, 2, 0, 0, 0, 0};

    // Check this with graph
    ASSERT_TRUE(adj_mat == adj_mat_truth);
}
