/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/io/frontend/utils/file_handle.hpp"
#include "detray/navigation/detail/helix.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/detail/whiteboard.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/types.hpp"
#include "detray/test/utils/detector_scanner.hpp"
#include "detray/test/utils/scan_utils.hpp"
#include "detray/test/utils/svg_display.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray::test {

/// @brief Test class that runs the helix scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_scan : public test::fixture_base<> {

    using algebra_t = typename detector_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using intersection_trace_t = typename detray::helix_scan<
        algebra_t>::template intersection_trace_type<detector_t>;
    // using track_generator_t =
    // uniform_track_generator<free_track_parameters_t>;
    using uniform_gen_t =
        random_numbers<scalar_t, std::uniform_real_distribution<scalar_t>>;
    using track_generator_t =
        random_track_generator<free_track_parameters_t, uniform_gen_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename track_generator_t::configuration;

        std::string m_name{"helix_scan"};
        // Save results for later use in downstream tests
        std::shared_ptr<test::whiteboard> m_white_board;
        bool m_write_inters{false};
        trk_gen_config_t m_trk_gen_cfg{};
        // Visualization style to be applied to the svgs
        detray::svgtools::styling::style m_style =
            detray::svgtools::styling::tableau_colorblind::style;

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
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
        config &whiteboard(std::shared_ptr<test::whiteboard> w_board) {
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
    explicit helix_scan(const detector_t &det,
                        const typename detector_t::name_map &names,
                        const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.whiteboard(cfg.whiteboard());
        m_cfg.write_intersections(cfg.write_intersections());
        m_cfg.track_generator() = cfg.track_generator();
    }

    /// Run the helix scan
    void TestBody() override {
        using nav_link_t = typename detector_t::surface_type::navigation_link;

        constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
        typename detector_t::geometry_context gctx{};

        // Index of the volume that the helix origin lies in
        dindex start_index{0u};

        // B-field vector for helix
        const typename fixture_type::point3 B{0.f * unit<scalar_t>::T,
                                              0.f * unit<scalar_t>::T,
                                              2.f * unit<scalar_t>::T};

        // Iterate through uniformly distributed momentum directions
        std::size_t n_tracks{0u};
        auto trk_state_generator = track_generator_t(m_cfg.track_generator());

        detray::io::file_handle outfile{
            m_cfg.name(), ".csv",
            std::ios::out | std::ios::binary | std::ios::trunc};

        if (m_cfg.write_intersections()) {
            *outfile << "index,type,x,y,z," << std::endl;
        }

        std::cout << "INFO: Running helix scan on: " << m_names.at(0) << "...\n"
                  << std::endl;

        for (auto trk : trk_state_generator) {

            // Get ground truth helix from track
            detail::helix helix(trk, &B);

            // Shoot helix through the detector and record all surfaces it
            // encounters
            const auto intersection_record =
                detector_scanner::run<detray::helix_scan>(
                    m_det, helix, 15.f * unit<scalar_t>::um, trk.p());

            // Csv output
            if (m_cfg.write_intersections()) {
                for (const auto &single_ir : intersection_record) {
                    const auto &intersection = single_ir.intersection;
                    const auto sf =
                        detray::surface{m_det, intersection.sf_desc};
                    auto glob_pos = sf.local_to_global(gctx, intersection.local,
                                                       helix.dir());
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
                detector_scanner::trace_intersections<leaving_world>(
                    intersection_record, start_index);

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

                detail::svg_display(gctx, il, intersection_record, helix,
                                    "helix", m_cfg.name());
            }

            // Is the succession of volumes consistent ?
            ASSERT_TRUE(err_code) << "\nFailed on helix " << n_tracks << "/"
                                  << trk_state_generator.size() << "\n"
                                  << helix;

            m_intersection_traces.push_back(std::move(intersection_record));

            ++n_tracks;
        }
        std::cout << "-----------------------------------\n"
                  << "Tested " << n_tracks << " helices: OK\n\n"
                  << "Adding " << m_intersection_traces.size()
                  << " helix traces to store\n"
                  << "-----------------------------------\n"
                  << std::endl;

        // Save the results
        m_cfg.whiteboard()->add("helix_scan", std::move(m_intersection_traces));
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
