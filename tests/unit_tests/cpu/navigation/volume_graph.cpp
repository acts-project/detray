/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project inlcude(s)
#include "detray/navigation/volume_graph.hpp"

#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/test/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <iostream>
#include <map>

// This tests the linking of a geometry by loading it into a graph structure
GTEST_TEST(detray_geometry, volume_graph) {
    using namespace detray;

    vecmem::host_memory_resource host_mr;

    toy_det_config toy_cfg{};
    toy_cfg.n_edc_layers(1u);

    auto [det, names] = create_toy_geometry(host_mr, toy_cfg);

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
