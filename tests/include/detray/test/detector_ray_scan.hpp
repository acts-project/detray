/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/utils/create_path.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/volume_graph.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/detail/whiteboard.hpp"
#include "detray/test/detector_scan_config.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/types.hpp"
#include "detray/test/utils/detector_scan_utils.hpp"
#include "detray/test/utils/detector_scanner.hpp"
#include "detray/test/utils/hash_tree.hpp"

// System include(s)
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace detray::test {

/// @brief Test class that runs the ray scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class ray_scan : public test::fixture_base<> {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using ray_t = detail::ray<algebra_t>;
    using intersection_trace_t = typename detray::ray_scan<
        algebra_t>::template intersection_trace_type<detector_t>;
    using uniform_gen_t =
        random_numbers<scalar_t, std::uniform_real_distribution<scalar_t>>;
    using track_generator_t = random_track_generator<ray_t, uniform_gen_t>;

    public:
    using fixture_type = test::fixture_base<>;
    using config = detector_scan_config<track_generator_t>;

    template <typename config_t>
    explicit ray_scan(const detector_t &det,
                      const typename detector_t::name_map &names,
                      const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.intersection_file(cfg.intersection_file());
        m_cfg.track_param_file(cfg.track_param_file());
        m_cfg.whiteboard(cfg.whiteboard());
        m_cfg.mask_tolerance(cfg.mask_tolerance());
        m_cfg.track_generator() = cfg.track_generator();
        m_cfg.write_intersections(cfg.write_intersections());
    }

    /// Run the ray scan
    void TestBody() override {

        // Get the volume adjaceny matrix from the detector
        volume_graph graph(m_det);
        const auto &adj_mat = graph.adjacency_matrix();

        // Fill adjacency matrix from ray scan and compare
        dvector<dindex> adj_mat_scan(adj_mat.size(), 0);
        // Keep track of the objects that have already been seen per volume
        std::unordered_set<dindex> obj_hashes = {};

        // Index of the volume that the ray origin lies in
        dindex start_index{0u};

        // Count the number of rays
        std::size_t n_tracks{0u};

        auto ray_generator = track_generator_t(m_cfg.track_generator());
        m_intersection_traces.reserve(ray_generator.size());

        std::cout << "\nINFO: Running ray scan on: " << m_names.at(0) << "\n"
                  << std::endl;

        if (io::file_exists(m_cfg.intersection_file()) &&
            io::file_exists(m_cfg.track_param_file())) {

            std::cout << "INFO: Reading data from file...\n" << std::endl;

            // Fill the intersection traces from file
            detector_scanner::read(m_cfg.intersection_file(),
                                   m_cfg.track_param_file(),
                                   m_intersection_traces);
        } else {

            std::cout << "INFO: Generating trace data...\n" << std::endl;

            for (const auto &ray : ray_generator) {

                // Record all intersections and surfaces along the ray
                auto intersection_trace =
                    detector_scanner::run<detray::ray_scan>(
                        m_det, ray, m_cfg.mask_tolerance());

                ASSERT_FALSE(intersection_trace.empty()) << ray;

                m_intersection_traces.push_back(std::move(intersection_trace));
            }
        }

        std::cout << "INFO: Checking trace data...\n" << std::endl;

        // Iterate through the scan data and perfrom checks
        for (const auto &intersection_trace : m_intersection_traces) {

            assert((intersection_trace.size() > 0) &&
                   "Invalid intersection trace");

            // Retrieve the test ray
            const auto &track = intersection_trace.front().track_param;
            detail::ray<algebra_t> ray(track);

            // Run consistency checks on the trace
            bool success = detector_scanner::check_trace<detector_t>(
                intersection_trace, start_index, adj_mat_scan, obj_hashes);

            // Display the detector, track and intersections for debugging
            if (not success) {
                detector_scanner::display_error(
                    m_gctx, m_det, m_names, m_cfg.name(), ray,
                    intersection_trace, m_cfg.svg_style(), n_tracks,
                    ray_generator.size());
            }

            ASSERT_TRUE(success);

            ++n_tracks;
        }
        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " rays: OK\n\n"
                  << "Adding " << m_intersection_traces.size()
                  << " ray traces to store\n";

        // Save the results

        // Csv output
        if (m_cfg.write_intersections()) {
            detector_scanner::write(m_cfg.name() + "_intersections.csv",
                                    m_cfg.name() + "_track_parameters.csv",
                                    m_intersection_traces);

            std::cout << "\nWrote  " << m_intersection_traces.size()
                      << " ray traces to file\n";
        }

        // Move the data to the whiteboard
        m_cfg.whiteboard()->add("ray_scan", std::move(m_intersection_traces));

        std::cout << "-----------------------------------\n" << std::endl;

        // Check that the links that were discovered by the scan match the
        // volume graph
        // ASSERT_TRUE(adj_mat == adj_mat_scan) <<
        // detector_scanner::print_adj(adj_mat_scan);

        // Compare the adjacency that was discovered in the ray scan to the
        // hashed one for the toy detector.
        // The hash tree is still Work in Progress !
        /*auto geo_checker = hash_tree(adj_mat);
        const bool check_links = (geo_checker.root() == root_hash);

        std::cout << "All links reachable: " << (check_links ? "OK" : "FAILURE")
                << std::endl;*/
    }

    private:
    /// The configuration of this test
    config m_cfg;
    /// The geometry context to scan
    typename detector_t::geometry_context m_gctx{};
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
    /// Record the scan for later use
    std::vector<intersection_trace_t> m_intersection_traces;
};

}  // namespace detray::test
