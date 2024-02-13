/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/surface.hpp"
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/io/frontend/utils/file_handle.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/types.hpp"
#include "detray/test/utils/particle_gun.hpp"
#include "detray/test/utils/ray_scan_utils.hpp"
#include "detray/test/utils/svg_display.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the helix scan on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_scan : public test::fixture_base<> {

    using transform3_t = typename detector_t::transform3;
    using free_track_parameters_t = free_track_parameters<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename uniform_track_generator<
            free_track_parameters_t>::configuration;

        std::string m_name{"helix_scan"};
        bool m_write_inters{false};
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
    explicit helix_scan(const detector_t &det,
                        const typename detector_t::name_map &names,
                        const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.write_intersections(cfg.write_intersections());
        m_cfg.track_generator() = cfg.track_generator();
    }

    /// Run the helix scan
    void TestBody() override {
        using scalar_t = typename detector_t::scalar_type;
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
        auto trk_state_generator =
            uniform_track_generator<free_track_parameters_t>(
                m_cfg.track_generator());

        detray::io::file_handle outfile{
            m_cfg.name(), ".csv",
            std::ios::out | std::ios::binary | std::ios::trunc};

        if (m_cfg.write_intersections()) {
            *outfile << "index,type,x,y,z," << std::endl;
        }

        std::cout << "INFO: Running helix scan on: " << m_names.at(0) << "\n("
                  << trk_state_generator.size() << " helices) ...\n"
                  << std::endl;

        for (auto trk : trk_state_generator) {

            // Get ground truth helix from track
            detail::helix helix(trk, &B);

            // Shoot helix through the detector and record all surfaces it
            // encounters
            const auto intersection_record = particle_gun::shoot_particle(
                m_det, helix, 14.9f * unit<scalar_t>::um);

            // Csv output
            if (m_cfg.write_intersections()) {
                for (const auto &single_ir : intersection_record) {
                    const auto &intersection = single_ir.second;
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

                detail::svg_display(gctx, il, intersection_record, helix,
                                    "helix", m_cfg.name());
            }

            // Is the succession of volumes consistent ?
            ASSERT_TRUE(err_code) << "\nFailed on helix " << n_tracks << "/"
                                  << trk_state_generator.size() << "\n"
                                  << helix;

            ++n_tracks;
        }
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
