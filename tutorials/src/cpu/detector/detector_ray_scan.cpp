/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/build_toy_detector.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/volume_graph.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/utils/hash_tree.hpp"
#include "detray/test/utils/particle_gun.hpp"
#include "detray/test/utils/ray_scan_utils.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>

// Hash of the "correct" geometry
constexpr std::size_t root_hash = 3244;

/// Check a given detecor for consistent linking by shooting rays/helices and
/// recording every intersection with the geometry. This intersection record
/// can then be checked for matching portals at the volume boundary surfaces
/// ( @c trace_intersections ) and checked for a consistent 'path' from volume
/// to volume ( @c check_consistency ). See also documentation in
/// 'tests/common/tools/ray_scan_utils.hpp'.
int main() {

    // Can also be performed with helices
    using ray_t = detray::detail::ray<detray::tutorial::algebra_t>;

    // Build the geometry
    vecmem::host_memory_resource host_mr;
    const auto [det, names] = detray::build_toy_detector(host_mr);

    // The invalid link value for the toy detector
    using nav_link_t = typename decltype(det)::surface_type::navigation_link;
    constexpr auto leaving_world{detray::detail::invalid_value<nav_link_t>()};

    // Get the volume adjaceny matrix from ray scan
    detray::volume_graph graph(det);
    const auto &adj_mat = graph.adjacency_matrix();  // < need this for the size
    detray::dvector<detray::dindex> adj_mat_scan(adj_mat.size(), 0);

    // Keep track of the objects that have already been seen per volume
    std::unordered_set<detray::dindex> obj_hashes = {};

    // Index of the volume that the ray origin lies in
    detray::dindex start_index{0u};

    // Generate a number of rays
    auto ray_generator = detray::uniform_track_generator<ray_t>{};
    ray_generator.config().theta_steps(100u).phi_steps(100u);

    // Run the check
    std::cout << "\nScanning detector (" << ray_generator.size()
              << " rays) ...\n"
              << std::endl;
    bool success = true;
    for (const auto ray : ray_generator) {

        // Record all intersections and surfaces along the ray
        const auto intersection_record =
            detray::particle_gun::shoot_particle(det, ray);

        // Create a trace of the volume indices that were encountered
        // and check that portal intersections are connected
        auto [portal_trace, surface_trace, err_code] =
            detray::trace_intersections<leaving_world>(intersection_record,
                                                       start_index);
        success &= err_code;

        // Is the succession of volumes consistent ?
        success &= detray::check_connectivity<leaving_world>(portal_trace);

        // Build an adjacency matrix from this trace that can be checked against
        // the geometry hash (see 'track_geometry_changes')
        detray::build_adjacency<leaving_world>(portal_trace, surface_trace,
                                               adj_mat_scan, obj_hashes);
    }

    // Check result
    std::cout << "Ray scan: " << (success ? "OK" : "FAILURE") << std::endl;

    // Compare the adjacency that was discovered in the ray scan to the hashed
    // one for the toy detector.
    // The hash tree is still Work in Progress !
    auto geo_checker = detray::hash_tree(adj_mat);
    const bool check_links = (geo_checker.root() == root_hash);

    std::cout << "All links reachable: " << (check_links ? "OK" : "FAILURE")
              << std::endl;
}
