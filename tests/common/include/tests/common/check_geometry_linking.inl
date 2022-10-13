/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/geometry/volume_graph.hpp"
#include "tests/common/tools/hash_tree.hpp"
//#include "tests/common/tools/read_geometry.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

// vecmem::host_memory_resource host_mr;
// auto [d, name_map] = read_from_csv(tml_files, host_mr);

namespace __plugin {

constexpr std::size_t root_hash = 3244;

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, check_geometry_linking) {

    // Detector configuration
    constexpr std::size_t n_brl_layers{4};
    constexpr std::size_t n_edc_layers{3};
    vecmem::host_memory_resource host_mr;
    auto det = create_toy_geometry(host_mr, n_brl_layers, n_edc_layers);

    // Build the graph
    volume_graph graph(det);
    const auto &adj_mat = graph.adjacency_matrix();

    // std::cout << graph.to_string() << std::endl;

    auto geo_checker = hash_tree(adj_mat);

    EXPECT_EQ(geo_checker.root(), root_hash);
}

}  // namespace __plugin
