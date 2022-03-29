/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include <iostream>
#include <vecmem/memory/host_memory_resource.hpp>

#include "detray/geometry/volume_graph.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/hash_tree.hpp"
#include "tests/common/tools/ray_scan_utils.hpp"

/// @note __plugin has to be defined with a preprocessor command
using namespace detray;

constexpr std::size_t root_hash =
    687;  // TODO: Find hash function wihtout coll.!

/** Print and adjacency list */
void print_adj(const dvector<dindex> &adjacency_matrix) {

    std::size_t dim = static_cast<dindex>(std::sqrt(adjacency_matrix.size()));

    for (std::size_t i = 0; i < dim - 1; ++i) {
        std::cout << "[>>] Node with index " << i << std::endl;
        std::cout << " -> edges: " << std::endl;
        for (std::size_t j = 0; j < dim; ++j) {
            const auto degr = adjacency_matrix[dim * i + j];
            if (degr == 0) {
                continue;
            }
            std::string n_occur =
                degr > 1 ? "\t\t\t\t(" + std::to_string(degr) + "x)" : "";

            // Edge that leads out of the detector world
            if (j == dim - 1 and degr != 0) {
                std::cout << "    -> leaving world " + n_occur << std::endl;
            } else {
                std::cout << "    -> " << std::to_string(j) + "\t" + n_occur
                          << std::endl;
            }
        }
    }
}

// Tests the consistency of the toy geometry implementation. In principle,
// every geometry can be checked this way.
TEST(ALGEBRA_PLUGIN, geometry_scan) {

    // Build the geometry (modeled as a unified index geometry)
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);

    // Build the graph
    volume_graph graph(toy_det);
    const auto &adj_mat = graph.adjacency_matrix();

    // std::cout << graph.to_string() << std::endl;

    // Now get the adjaceny matrix from ray scan
    dvector<dindex> adj_mat_scan(adj_mat.size(), 0);
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    unsigned int theta_steps = 100;
    unsigned int phi_steps = 100;
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
            build_adjacency(portal_trace, surface_trace, adj_mat_scan,
                            obj_hashes);
        }
    }

    // print_adj(adj_scan);
    // ASSERT_TRUE(adj_mat == adj_mat_truth);

    // ASSERT_EQ(adj_linking, adj_scan);
    // auto geo_checker = hash_tree(adj_mat_scan);
    // ASSERT_EQ(geo_checker.root(), root_hash);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
