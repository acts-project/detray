/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/tools/geometry_graph.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/hash_tree.hpp"
#include "tests/common/tools/ray_scan_utils.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

constexpr std::size_t vol0_hash = 2;
constexpr std::size_t vol1_hash = 2;  // TODO: Find hash function wihtout coll.!
constexpr std::size_t vol2_hash = 10;
constexpr std::size_t vol3_hash = 1798;

/** Print and adjacency list */
void print_adj(std::map<dindex, std::map<dindex, dindex>> &adjacency_list) {
    const auto print_neighbor =
        [&](const std::pair<const dindex, const dindex> &n) -> std::string {
        // Print the number of edges, if it is greater than one
        std::string n_occur =
            n.second > 1 ? "\t\t(" + std::to_string(n.second) + "x)" : "";

        // Edge that leads out of the detector world
        if (n.first == dindex_invalid) {
            return "leaving world" + n_occur;
        }

        return std::to_string(n.first) + "\t\t\t" + n_occur;
    };

    for (const auto &[vol_index, neighbors] : adjacency_list) {
        std::cout << "[>>] Node with index " << vol_index << std::endl;
        std::cout << " -> neighbors: " << std::endl;
        for (const auto &nbr : neighbors) {
            std::cout << "    -> " << print_neighbor(nbr) << std::endl;
        }
    }
}

// Tests the consistency of the toy geometry implementation. In principle,
// every geometry can be checked this way.
TEST(ALGEBRA_PLUGIN, geometry_scan) {

    // Build the geometry (modeled as a unified index geometry)
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);

    // Now get the adjaceny list from ray scan

    // Adjacency list to be filled in ray scan
    std::map<dindex, std::map<dindex, dindex>> adj_scan = {};
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    unsigned int theta_steps = 1000;
    unsigned int phi_steps = 1000;
    const point3 ori{0., 0., 0.};
    dindex start_index = 0;

    // Loops of theta values ]0,pi[
    for (unsigned int itheta = 0; itheta < theta_steps; ++itheta) {
        scalar theta = 0.001 + itheta * (M_PI - 0.001) / theta_steps;
        scalar sin_theta = std::sin(theta);
        scalar cos_theta = std::cos(theta);

        // Loops of phi values [-pi, pi]
        for (unsigned int iphi = 0; iphi < phi_steps; ++iphi) {
            // The direction
            scalar phi = -M_PI + iphi * (2 * M_PI) / phi_steps;
            scalar sin_phi = std::sin(phi);
            scalar cos_phi = std::cos(phi);
            const point3 dir{cos_phi * sin_theta, sin_phi * sin_theta,
                             cos_theta};

            // Record all intersections and objects along the ray
            const auto intersection_record = shoot_ray(toy_det, ori, dir);

            // Create a trace of the volume indices that were encountered
            auto [portal_trace, surface_trace] =
                trace_intersections(intersection_record, start_index);

            // Is this a sensible trace to be further examined?
            ASSERT_TRUE(check_connectivity(portal_trace));

            // Discover new linking information from this trace
            build_adjacency(portal_trace, surface_trace, adj_scan, obj_hashes);
        }
    }

    print_adj(adj_scan);

    // TODO: Join these sub trees into a single comprehensive tree
    /*auto geo_checker_vol0 =
        hash_tree<decltype(adj_scan.at(0)), dindex>(adj_scan.at(0));

    EXPECT_EQ(geo_checker_vol0.root(), vol0_hash);

    // This one fails, because the ray scan is kept very coarse for performance
    // reasons (run on the CI)
    auto geo_checker_vol1 =
        hash_tree<decltype(adj_scan.at(1)), dindex>(adj_scan.at(1));

    EXPECT_EQ(geo_checker_vol1.root(), vol1_hash);

    auto geo_checker_vol2 =
        hash_tree<decltype(adj_scan.at(2)), dindex>(adj_scan.at(2));

    EXPECT_EQ(geo_checker_vol2.root(), vol2_hash);

    auto geo_checker_vol3 =
        hash_tree<decltype(adj_scan.at(3)), dindex>(adj_scan.at(3));

    EXPECT_EQ(geo_checker_vol3.root(), vol3_hash);*/
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
