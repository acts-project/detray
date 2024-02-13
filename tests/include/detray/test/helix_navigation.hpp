/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/detectors/bfield.hpp"
#include "detray/navigation/detail/helix.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/test/fixture_base.hpp"
#include "detray/test/utils/navigation_check_helper.hpp"
#include "detray/test/utils/particle_gun.hpp"
#include "detray/test/utils/svg_display.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the helix navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_navigation : public test::fixture_base<> {

    using scalar_t = typename detector_t::scalar_type;
    using transform3_t = typename detector_t::transform3;
    using free_track_parameters_t = free_track_parameters<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename uniform_track_generator<
            free_track_parameters_t>::configuration;

        std::string m_name{"helix_navigation"};
        trk_gen_config_t m_trk_gen_cfg{};
        // Visualization style to be applied to the svgs
        detray::svgtools::styling::style m_style =
            detray::svgtools::styling::tableau_colorblind::style;

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
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
        /// @}
    };

    template <typename config_t>
    explicit helix_navigation(const detector_t &det,
                              const typename detector_t::name_map &names,
                              const config_t &cfg = {})
        : m_det{det}, m_names{names} {
        m_cfg.name(cfg.name());
        m_cfg.track_generator() = cfg.track_generator();
        m_cfg.propagation() = cfg.propagation();
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        /// Type that holds the intersection information
        using intersection_t =
            intersection2D<typename detector_t::surface_type, transform3_t>;

        /// Inspector that records all encountered surfaces
        using object_tracer_t =
            navigation::object_tracer<intersection_t, dvector,
                                      navigation::status::e_on_module,
                                      navigation::status::e_on_portal>;
        /// Inspector that prints the navigator state from within the
        /// navigator's method calls (cannot be done with an actor)
        using nav_print_inspector_t = navigation::print_inspector;
        /// Aggregation of multiple inspectors
        using inspector_t =
            aggregate_inspector<object_tracer_t, nav_print_inspector_t>;

        // Navigation with inspection
        using navigator_t = navigator<detector_t, inspector_t, intersection_t>;
        // Runge-Kutta stepper
        using bfield_t = bfield::const_field_t;
        using stepper_t = rk_stepper<typename bfield_t::view_t, transform3_t,
                                     unconstrained_step, stepper_rk_policy,
                                     stepping::print_inspector>;
        // Propagator with pathlimit aborter
        using actor_chain_t = actor_chain<dtuple, pathlimit_aborter>;
        using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

        // Use default context
        typename detector_t::geometry_context gctx{};

        // Propagator
        propagator_t prop{m_cfg.propagation()};

        // B-field vector for helix
        const typename fixture_type::point3 B{0.f * unit<scalar_t>::T,
                                              0.f * unit<scalar_t>::T,
                                              2.f * unit<scalar_t>::T};
        bfield_t hom_bfield = bfield::create_const_field(B);

        // Iterate through uniformly distributed momentum directions
        std::size_t n_tracks{0u};
        auto trk_state_generator =
            uniform_track_generator<free_track_parameters_t>(
                m_cfg.track_generator());

        std::cout << "INFO: Running helix navigation check on: "
                  << m_names.at(0) << "\n(" << trk_state_generator.size()
                  << " helices) ...\n"
                  << std::endl;

        /// Error statistic
        std::size_t n_close_miss{0u}, n_fatal{0u};

        std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
        detray::io::file_handle debug_file{"./helix_navigation.txt", io_mode};

        for (auto track : trk_state_generator) {

            // Get ground truth helix from track
            detail::helix helix(track, &B);

            // Shoot helix through the detector and record all surfaces it
            // encounters
            const auto intersection_trace = particle_gun::shoot_particle(
                m_det, helix, 15.f * unit<scalar_t>::um);

            // Build actor and propagator states
            pathlimit_aborter::state pathlimit_aborter_state{
                m_cfg.propagation().stepping.path_limit};
            auto actor_states = std::tie(pathlimit_aborter_state);

            typename propagator_t::state propagation(track, hom_bfield, m_det);

            // Access to navigation information
            auto &nav_inspector = propagation._navigation.inspector();
            auto &obj_tracer = nav_inspector.template get<object_tracer_t>();
            auto &nav_printer = nav_inspector.template get<print_inspector>();

            // Acces to the stepper information
            auto &step_printer = propagation._stepping.inspector();

            // Run the propagation
            bool success = prop.propagate(propagation, actor_states);

            if (success) {
                success &= detail::compare_traces(intersection_trace,
                                                  obj_tracer, helix, n_tracks,
                                                  trk_state_generator.size());
            }
            if (not success) {
                // Write debug info to file
                *debug_file << "HELIX " << n_tracks << ":\n\n"
                            << nav_printer.to_string()
                            << step_printer.to_string();

                // Create the svg for failed tracks.
                detray::svgtools::illustrator il{m_det, m_names,
                                                 m_cfg.svg_style()};
                il.show_info(true);
                il.hide_eta_lines(true);
                il.hide_portals(false);
                il.hide_passives(false);

                detail::svg_display(gctx, il, intersection_trace, helix,
                                    "helix_" + std::to_string(n_tracks),
                                    m_cfg.name(), obj_tracer.object_trace);

                // Keep a statistic on the errors that occured
                if (!propagation._navigation.is_complete()) {
                    ++n_fatal;
                } else {
                    // @TODO: Check mask boundaries
                    ++n_close_miss;
                }
            }

            EXPECT_TRUE(success)
                << "\nFailed on helix " << n_tracks << "/"
                << trk_state_generator.size() << ": " << helix << "\n\n";

            ++n_tracks;
        }

        if (n_close_miss > 0u || n_fatal > 0u) {
            std::cout << "-----------------------------------"
                      << "Error Statistic:\n\n"
                      << "\n (close misses: " << n_close_miss
                      << ", fatal failures: " << n_fatal << ")\n"
                      << "-----------------------------------\n"
                      << std::endl;
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
