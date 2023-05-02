/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <map>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/volume_graph.hpp"

/// @note __plugin has to be defined with a preprocessor command

// This tests the linking of a geometry by loading it into a graph structure
TEST(ALGEBRA_PLUGIN, volume_graph) {
    using namespace detray;
    using namespace __plugin;

    vecmem::host_memory_resource host_mr;

    unsigned int n_brl_layers{4u};
    unsigned int n_edc_layers{1u};

    auto det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);
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

    // std::cout << graph.to_string() << std::endl;
    // std::cout << "Walking through geometry: " << std::endl;
    //  graph.bfs();

    const auto &adj_mat = graph.adjacency_matrix();

    // toy geometry 1 endcap layer, 4 barrel layers, including gap layers
    dvector<dindex> adj_mat_truth = {// beampipe
                                     1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 2,
                                     // endcap layer 1
                                     1, 108, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
                                     // connector layer
                                     1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1,
                                     // barrel layer 1
                                     1, 0, 1, 224, 1, 0, 0, 0, 0, 0, 1, 0, 0,
                                     // barrel gap 1
                                     0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0,
                                     // barrel layer 2
                                     0, 0, 1, 0, 1, 448, 1, 0, 0, 0, 1, 0, 0,
                                     // barrel gap 2
                                     0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0,
                                     // barrel layer 3
                                     0, 0, 1, 0, 0, 0, 1, 728, 1, 0, 1, 0, 0,
                                     // barrel gap 3
                                     0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0,
                                     // barrel layer 4
                                     0, 0, 1, 0, 0, 0, 0, 0, 1, 1092, 1, 0, 1,
                                     // connector layer
                                     1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
                                     // endcap laer 2
                                     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 108, 2,
                                     // world volume
                                     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // Check this with graph
    ASSERT_TRUE(adj_mat == adj_mat_truth);
}
