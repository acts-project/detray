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
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/detail/whiteboard.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/types.hpp"
#include "detray/test/utils/detector_scanner.hpp"
#include "detray/test/utils/hash_tree.hpp"
#include "detray/test/utils/scan_utils.hpp"
#include "detray/test/utils/svg_display.hpp"

// System include(s)
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace detray::test {

/// Print and adjacency list
inline std::string print_adj(const dvector<dindex> &adjacency_matrix) {

    std::size_t dim = static_cast<dindex>(math::sqrt(adjacency_matrix.size()));
    std::stringstream out_stream{};

    for (std::size_t i = 0u; i < dim - 1; ++i) {
        out_stream << "[>>] Node with index " << i << std::endl;
        out_stream << " -> edges: " << std::endl;
        for (std::size_t j = 0u; j < dim; ++j) {
            const auto degr = adjacency_matrix[dim * i + j];
            if (degr == 0) {
                continue;
            }
            std::string n_occur =
                degr > 1 ? "\t\t\t\t(" + std::to_string(degr) + "x)" : "";

            // Edge that leads out of the detector world
            if (j == dim - 1 and degr != 0) {
                out_stream << "    -> leaving world " + n_occur << std::endl;
            } else {
                out_stream << "    -> " << std::to_string(j) + "\t" + n_occur
                           << std::endl;
            }
        }
    }

    return out_stream.str();
}

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

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename track_generator_t::configuration;

        std::string m_name{"ray_scan"};
        // Save results for later use in downstream tests
        std::shared_ptr<test::whiteboard> m_white_board;
        // Name of the input file, containing the complete ray scan traces
        std::string m_intersection_file{""};
        std::string m_track_param_file{""};
        // Mask tolerance for the intersectors
        std::array<scalar_t, 2> m_mask_tol{
            std::numeric_limits<scalar_t>::epsilon(),
            std::numeric_limits<scalar_t>::epsilon()};
        // Configuration of the ray generator
        trk_gen_config_t m_trk_gen_cfg{};
        // Write intersection points for plotting
        bool m_write_inters{false};
        // Visualization style to be applied to the svgs
        detray::svgtools::styling::style m_style =
            detray::svgtools::styling::tableau_colorblind::style;

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        const std::string &intersection_file() const {
            return m_intersection_file;
        }
        const std::string &track_param_file() const {
            return m_track_param_file;
        }
        std::array<scalar_t, 2> mask_tolerance() const { return m_mask_tol; }
        std::shared_ptr<test::whiteboard> whiteboard() { return m_white_board; }
        std::shared_ptr<test::whiteboard> whiteboard() const {
            return m_white_board;
        }
        bool write_intersections() const { return m_write_inters; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
        const auto &svg_style() const { return m_style; }
        /// @}

        /// Setters
        /// @{
        config &name(const std::string n) {
            m_name = n;
            return *this;
        }
        config &intersection_file(const std::string f) {
            m_intersection_file = std::move(f);
            return *this;
        }
        config &track_param_file(const std::string f) {
            m_track_param_file = std::move(f);
            return *this;
        }
        config &mask_tolerance(const std::array<scalar_t, 2> tol) {
            m_mask_tol = tol;
            return *this;
        }
        config &whiteboard(std::shared_ptr<test::whiteboard> w_board) {
            if (!w_board) {
                throw std::invalid_argument(
                    "Ray scan: No valid whiteboard instance");
            }
            m_white_board = std::move(w_board);
            return *this;
        }
        config &write_intersections(const bool do_write) {
            m_write_inters = do_write;
            return *this;
        }
        /// @}
    };

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
        using nav_link_t = typename detector_t::surface_type::navigation_link;

        constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
        typename detector_t::geometry_context gctx{};

        // Get the volume adjaceny matrix from ray scan
        volume_graph graph(m_det);
        const auto &adj_mat = graph.adjacency_matrix();
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
            auto track(intersection_trace.front().track_param);
            detail::ray<algebra_t> ray(track);

            // Create a trace of the volume indices that were encountered
            // and check that portal intersections are connected
            auto [portal_trace, surface_trace, err_code] =
                detector_scanner::trace_intersections<leaving_world>(
                    intersection_trace, start_index);

            // Is the succession of volumes consistent ?
            err_code &= detector_scanner::check_connectivity<leaving_world>(
                portal_trace);

            // Display the detector, track and intersections for debugging
            if (not err_code) {

                // Creating the svg generator for the detector.
                detray::svgtools::illustrator il{m_det, m_names,
                                                 m_cfg.svg_style()};
                il.show_info(true);
                il.hide_eta_lines(true);
                il.hide_portals(false);
                il.hide_passives(false);

                detail::svg_display(gctx, il, intersection_trace, ray, "ray",
                                    m_cfg.name());
            }

            ASSERT_TRUE(err_code) << "\nFailed on ray " << n_tracks << "/"
                                  << ray_generator.size() << "\n"
                                  << ray;

            // Build an adjacency matrix from this trace that can be checked
            // against the geometry hash (see 'track_geometry_changes')
            detector_scanner::build_adjacency<leaving_world>(
                portal_trace, surface_trace, adj_mat_scan, obj_hashes);

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
        // ASSERT_TRUE(adj_mat == adj_mat_scan) << print_adj(adj_mat_scan);

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
    /// The detector to be checked
    const detector_t &m_det;
    /// Volume names
    const typename detector_t::name_map &m_names;
    /// Record the scan for later use
    std::vector<intersection_trace_t> m_intersection_traces;
};

}  // namespace detray::test
