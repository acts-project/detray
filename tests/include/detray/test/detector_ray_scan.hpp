/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/frontend/utils/file_handle.hpp"
#include "detray/navigation/detail/ray.hpp"
#include "detray/navigation/volume_graph.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/types.hpp"
#include "detray/test/utils/hash_tree.hpp"
#include "detray/test/utils/particle_gun.hpp"
#include "detray/test/utils/ray_scan_utils.hpp"
#include "detray/test/utils/svg_display.hpp"

// System include(s)
#include <iostream>
#include <sstream>
#include <string>

namespace detray {

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

    using transform3_t = typename detector_t::transform3;
    using ray_t = detail::ray<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t =
            typename uniform_track_generator<ray_t>::configuration;

        std::string m_name{"ray_scan"};
        // Write intersection points for plotting
        bool m_write_inters{false};
        // Configuration of the ray generator
        trk_gen_config_t m_trk_gen_cfg{};
        // Visualization style to be applied to the svgs
        detray::svgtools::styling::style m_style =
            detray::svgtools::styling::tableau_colorblind::style;

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
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
        m_cfg.write_intersections(cfg.write_intersections());
        m_cfg.track_generator() = cfg.track_generator();
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

        std::size_t n_tracks{0u};
        auto ray_generator =
            uniform_track_generator<ray_t>(m_cfg.track_generator());

        // Csv output file
        detray::io::file_handle outfile{
            m_cfg.name(), ".csv",
            std::ios::out | std::ios::binary | std::ios::trunc};

        if (m_cfg.write_intersections()) {
            *outfile << "index,type,x,y,z," << std::endl;
        }

        std::cout << "INFO: Running ray scan on: " << m_names.at(0) << "\n("
                  << ray_generator.size() << " rays) ...\n"
                  << std::endl;

        for (const auto &ray : ray_generator) {

            // Record all intersections and surfaces along the ray
            const auto intersection_record =
                particle_gun::shoot_particle(m_det, ray);

            // Csv output
            if (m_cfg.write_intersections()) {
                for (const auto &single_ir : intersection_record) {
                    const auto &intersection = single_ir.second;
                    const auto sf =
                        detray::surface{m_det, intersection.sf_desc};
                    auto glob_pos =
                        sf.local_to_global(gctx, intersection.local, ray.dir());
                    *outfile
                        << n_tracks << ","
                        << static_cast<int>(intersection.sf_desc.barcode().id())
                        << "," << glob_pos[0] << "," << glob_pos[1] << ","
                        << glob_pos[2] << "," << std::endl;
                }
            }

            // Create a trace of the volume indices that were encountered
            // and check that portal intersections are connected
            auto [portal_trace, surface_trace, err_code] =
                trace_intersections<leaving_world>(intersection_record,
                                                   start_index);

            // Is the succession of volumes consistent ?
            err_code &= check_connectivity<leaving_world>(portal_trace);

            // Display the detector, track and intersections for debugging
            if (not err_code) {

                // Creating the svg generator for the detector.
                detray::svgtools::illustrator il{m_det, m_names,
                                                 m_cfg.svg_style()};
                il.show_info(true);
                il.hide_eta_lines(true);
                il.hide_portals(false);
                il.hide_passives(false);

                detail::svg_display(gctx, il, intersection_record, ray, "ray",
                                    m_cfg.name());
            }

            ASSERT_TRUE(err_code) << "\nFailed on ray " << n_tracks << "/"
                                  << ray_generator.size() << "\n"
                                  << ray;

            // Build an adjacency matrix from this trace that can be checked
            // against the geometry hash (see 'track_geometry_changes')
            build_adjacency<leaving_world>(portal_trace, surface_trace,
                                           adj_mat_scan, obj_hashes);

            ++n_tracks;
        }

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
};

}  // namespace detray
