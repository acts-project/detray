/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/geometry/volume_graph.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/io/json/json_reader.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "tests/common/tools/hash_tree.hpp"
#include "tests/common/tools/particle_gun.hpp"
#include "tests/common/tools/ray_scan_utils.hpp"

// Example linear algebra plugin: std::array
#include "detray/tutorial/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <fstream>
#include <iostream>
#include <string>

// Hash of the "correct" geometry
constexpr std::size_t root_hash = 3244;

/// Execution function for the ray scan
///
/// @tparam detector_t the type of the detector (toy, read in)
/// @param det the detector object
/// @param theta_steps number of steps in theta
/// @param phi_steps number of steps in phi
///
/// @return 1 if the scan was successful, 0 otherwise
template <typename detector_t>
int run_ray_scan(const detector_t &det, unsigned int theta_steps = 100u,
                 unsigned int phi_steps = 100u,
                 const std::string &outfile = "") {

    std::fstream out;
    if (not outfile.empty()) {
        out.open(outfile, std::ios::out);
        out << "ray_index,x,y,z," << std::endl;
    }

    // Can also be performed with helices
    using ray_t = detray::detail::ray<detray::tutorial::transform3>;

    std::cout << "Number of rays in theta: " << theta_steps << std::endl;
    std::cout << "Number of rays in phi: " << phi_steps << std::endl;

    // The invalid link value for the toy detector
    using nav_link_t = typename detector_t::surface_type::navigation_link;
    constexpr auto leaving_world{detray::detail::invalid_value<nav_link_t>()};

    // Get the volume adjaceny matrix from ray scan
    detray::volume_graph graph(det);
    const auto &adj_mat = graph.adjacency_matrix();  // < need this for the size
    detray::dvector<detray::dindex> adj_mat_scan(adj_mat.size(), 0);

    // Keep track of the objects that have already been seen per volume
    std::unordered_set<detray::dindex> obj_hashes = {};

    // Index of the volume that the ray origin lies in
    detray::dindex start_index{0u};

    // Origin of the rays
    const detray::tutorial::point3 origin{0.f, 0.f, 0.f};
    auto ray_generator =
        detray::uniform_track_generator<ray_t>(theta_steps, phi_steps, origin);

    // Geometry context for the ray scan & ray counter
    typename detector_t::geometry_context gctx{};
    std::size_t ray_counter{0u};

    // Run the check
    std::cout << "\nScanning detector (" << ray_generator.size()
              << " rays) ...\n"
              << std::endl;
    bool success = true;
    for (const auto ray : ray_generator) {

        // Record all intersections and surfaces along the ray
        const auto intersection_record =
            detray::particle_gun::shoot_particle(det, ray);

        // Write the intersection record to file if requested
        if (not outfile.empty()) {
            for (const auto &single_ir : intersection_record) {
                const auto &intersection = single_ir.second;
                const auto &sf = detray::surface{det, intersection.surface};
                auto glob_pos = sf.local_to_global(gctx, intersection.local);
                out << ray_counter << "," << glob_pos[0] << "," << glob_pos[1]
                    << "," << glob_pos[2] << "," << std::endl;
            }
        }

        // Create a trace of the volume indices that were encountered
        // and check that portal intersections are connected
        auto [portal_trace, surface_trace] =
            detray::trace_intersections<leaving_world>(intersection_record,
                                                       start_index);

        // Is the succession of volumes consistent ?
        success &= detray::check_connectivity<leaving_world>(portal_trace);

        // Build an adjacency matrix from this trace that can be checked against
        // the geometry hash (see 'track_geometry_changes')
        detray::build_adjacency<leaving_world>(portal_trace, surface_trace,
                                               adj_mat_scan, obj_hashes);

         // Increase the ray counter
         ++ray_counter;                                      
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
    if (not outfile.empty()) {
        out.close();
    }
    return 1;
}

/// Check a given detecor for consistent linking by shooting rays/helices and
/// recording every intersection with the geometry. This intersection record
/// can then be checked for matching portals at the volume boundary surfaces
/// ( @c trace_intersections ) and checked for a consistent 'path' from volume
/// to volume ( @c check_consistency ). See also documentation in
/// 'tests/common/tools/ray_scan_utils.hpp'.
/// Check a given detecor for consistent linking by shooting rays/helices and
/// recording every intersection with the geometry. This intersection record
/// can then be checked for matching portals at the volume boundary surfaces
/// ( @c trace_intersections ) and checked for a consistent 'path' from volume
/// to volume ( @c check_consistency ). See also documentation in
/// 'tests/common/tools/ray_scan_utils.hpp'.
int main(int argc, char *argv[]) {

    vecmem::host_memory_resource host_mr;

    // Input steering of the main program, file name of the detector
    if (argc < 2) {
        std::cout << "Running ray scan on ad-hoc built toy detector."
                  << std::endl;
        // Build the geometry
        auto det = detray::create_toy_geometry(host_mr);
        return run_ray_scan(det);
    } else if (argc < 4) {
        std::cout << "Running ray scan from command line arguments, requires"
                     " detector file name, theta and phi steps"
                  << std::endl;
        return -1;
    }

    // Command line arguments conversion
    std::string filename(argv[1]);
    auto theta_steps = static_cast<unsigned int>(std::stoi(argv[2]));
    auto phi_steps = static_cast<unsigned int>(std::stoi(argv[3]));
    std::string outfile = "";
    if (argc > 4) {
        outfile = argv[4];
    }

    // Read the detector in
    std::cout << "Running ray scan from input detector, reading file: "
              << argv[1] << std::endl;

    using detector_t = detray::detector<detray::default_metadata>;

    // @todo : Create volume name map in 'create_toy_detector'
    typename detector_t::name_map volume_name_map = {{0u, "json_detector"}};

    // Read-in detector
    detector_t det{host_mr};

    detray::json_geometry_reader<detector_t> geo_reader;
    geo_reader.read(det, volume_name_map, filename);

    return run_ray_scan(det, theta_steps, phi_steps, outfile);
}
