/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/types.hpp"
#include "tests/common/tools/hash_tree.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/ray_scan_utils.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Gtest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <iostream>

using namespace detray;

// constexpr std::size_t root_hash =
//     687u;  // TODO: Find hash function wihtout coll.!

using transform3_type = test::transform3;
using ray_type = detail::ray<transform3_type>;

/// Print and adjacency list
void print_adj(const dvector<dindex> &adjacency_matrix) {

    std::size_t dim = static_cast<dindex>(std::sqrt(adjacency_matrix.size()));

    for (std::size_t i = 0u; i < dim - 1; ++i) {
        std::cout << "[>>] Node with index " << i << std::endl;
        std::cout << " -> edges: " << std::endl;
        for (std::size_t j = 0u; j < dim; ++j) {
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
GTEST_TEST(detray_geometry, toy_geometry_scan) {

    // Build the geometry (modeled as a unified index geometry)
    vecmem::host_memory_resource host_mr;
    auto toy_det = create_toy_geometry(host_mr);
    using nav_link_t =
        typename decltype(toy_det)::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Build the graph
    volume_graph graph(toy_det);
    const auto &adj_mat = graph.adjacency_matrix();

    // std::cout << graph.to_string() << std::endl;

    // Now get the adjaceny matrix from ray scan
    dvector<dindex> adj_mat_scan(adj_mat.size(), 0);
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    unsigned int theta_steps{100u};
    unsigned int phi_steps{100u};
    const point3 ori{0.f, 0.f, 0.f};
    dindex start_index{0u};

    // Iterate through uniformly distributed momentum directions
    for (const auto test_ray :
         uniform_track_generator<ray_type>(theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            particle_gun::shoot_particle(toy_det, test_ray);

        // Create a trace of the volume indices that were encountered
        auto [portal_trace, surface_trace] = trace_intersections<leaving_world>(
            intersection_record, start_index);

        // Is this a sensible trace to be further examined?
        ASSERT_TRUE(check_connectivity<leaving_world>(portal_trace));

        // Discover new linking information from this trace
        build_adjacency<leaving_world>(portal_trace, surface_trace,
                                       adj_mat_scan, obj_hashes);
    }

    // print_adj(adj_scan);
    // ASSERT_TRUE(adj_mat == adj_mat_scan);

    // ASSERT_EQ(adj_linking, adj_scan);
    // auto geo_checker = hash_tree(adj_mat_scan);
    // ASSERT_EQ(geo_checker.root(), root_hash);
}

// @TODO: Create common check function for all detectors
/// Tests the consistency of the telescope geometry implementation.
GTEST_TEST(detray_geometry, telescope_geometry_scan) {

    vecmem::host_memory_resource host_mr;

    // Build telescope detector with 10 rectangles
    tel_det_config<rectangle2D<>> tel_cfg{20.f * unit<scalar>::mm,
                                          20.f * unit<scalar>::mm};
    tel_cfg.n_surfaces(10u).length(500.f * unit<scalar>::mm);

    const auto [tel_det, tel_names] =
        create_telescope_detector(host_mr, tel_cfg);

    using nav_link_t =
        typename decltype(tel_det)::surface_type::navigation_link;
    constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};

    // Build the graph
    volume_graph graph(tel_det);
    const auto &adj_mat = graph.adjacency_matrix();

    // For debugging:
    // std::cout << graph.to_string() << std::endl;

    // Now get the adjaceny matrix from ray scan
    dvector<dindex> adj_mat_scan(adj_mat.size(), 0);
    // Keep track of the objects that have already been seen per volume
    std::unordered_set<dindex> obj_hashes = {};

    unsigned int theta_steps{100u};
    unsigned int phi_steps{100u};
    const point3 ori{0.f, 0.f, -0.05f};
    dindex start_index{0u};

    // Iterate through uniformly distributed momentum directions
    for (const auto test_ray :
         uniform_track_generator<ray_type>(theta_steps, phi_steps, ori)) {

        // Record all intersections and objects along the ray
        const auto intersection_record =
            particle_gun::shoot_particle(tel_det, test_ray);

        // Create a trace of the volume indices that were encountered
        auto [portal_trace, surface_trace] = trace_intersections<leaving_world>(
            intersection_record, start_index);

        // Is this a sensible trace to be further examined?
        ASSERT_TRUE(check_connectivity<leaving_world>(portal_trace));

        // Discover new linking information from this trace
        build_adjacency<leaving_world>(portal_trace, surface_trace,
                                       adj_mat_scan, obj_hashes);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
