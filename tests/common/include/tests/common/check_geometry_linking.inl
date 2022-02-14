/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/tools/volume_graph.hpp"
#include "tests/common/tools/hash_tree.hpp"
//#include "tests/common/tools/read_geometry.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"

using namespace detray;

// vecmem::host_memory_resource host_mr;
// auto [d, name_map] = read_from_csv(tml_files, host_mr);

namespace __plugin {

constexpr std::size_t vol0_hash = 2;
constexpr std::size_t vol1_hash = 2;  // TODO: Find hash function without coll.!
constexpr std::size_t vol2_hash = 10;
constexpr std::size_t vol3_hash = 1798;

// This test runs intersection with all portals of the TrackML detector
TEST(ALGEBRA_PLUGIN, check_geometry_linking) {

    // Build the geometry (modeled as a unified index geometry)
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);
    using detector_t = decltype(toy_det);

    // Build the graph
    const auto g = volume_graph<detector_t>(
        toy_det.volumes(), toy_det.surfaces(), toy_det.mask_store());
    // const auto &adj_linking = g.adjacency_list();

    std::cout << g.to_string() << std::endl;

    // TODO: Join these sub trees into a single comprehensive tree
    /*auto geo_checker_vol0 =
        hash_tree<decltype(adj_linking.at(0)), dindex>(adj_linking.at(0));

    EXPECT_EQ(geo_checker_vol0.root(), vol0_hash);

    auto geo_checker_vol1 =
        hash_tree<decltype(adj_linking.at(1)), dindex>(adj_linking.at(1));

    EXPECT_EQ(geo_checker_vol1.root(), vol1_hash);

    auto geo_checker_vol2 =
        hash_tree<decltype(adj_linking.at(2)), dindex>(adj_linking.at(2));

    EXPECT_EQ(geo_checker_vol2.root(), vol2_hash);

    auto geo_checker_vol3 =
        hash_tree<decltype(adj_linking.at(3)), dindex>(adj_linking.at(3));

    EXPECT_EQ(geo_checker_vol3.root(), vol3_hash);*/
}

}  // namespace __plugin
