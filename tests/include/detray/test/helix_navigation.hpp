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
#include "detray/test/utils/detector_scanner.hpp"
#include "detray/test/utils/navigation_check_helper.hpp"
#include "detray/test/utils/svg_display.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray::test {

/// @brief Test class that runs the helix navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_navigation : public test::fixture_base<> {

    using scalar_t = typename detector_t::scalar_type;
    using algebra_t = typename detector_t::algebra_type;
    using free_track_parameters_t = free_track_parameters<algebra_t>;
    using intersection_trace_t = typename detray::helix_scan<
        algebra_t>::template intersection_trace_type<detector_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {

        std::string m_name{"helix_navigation"};
        // Access to truth data
        std::shared_ptr<test::whiteboard> m_white_board;
        // The number of test tracks to run
        std::size_t m_n_tracks{detray::detail::invalid_value<std::size_t>()};
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
        std::size_t n_tracks() const { return m_n_tracks; }
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
        config &n_tracks(std::size_t n) {
            m_n_tracks = n;
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
        m_cfg.whiteboard(cfg.whiteboard());
        m_cfg.n_tracks(cfg.n_tracks());
        m_cfg.propagation() = cfg.propagation();

        if (!m_cfg.whiteboard()) {
            throw std::invalid_argument(
                "No white board was passed to helix navigation test");
        }
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        /// Type that holds the intersection information
        using intersection_t =
            intersection2D<typename detector_t::surface_type, algebra_t>;

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
        using stepper_t =
            rk_stepper<typename bfield_t::view_t, algebra_t, unconstrained_step,
                       stepper_rk_policy, stepping::print_inspector>;
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
        if (!m_cfg.whiteboard()->exists("helix_scan")) {
            throw std::runtime_error(
                "White board is empty! Please run helix scan first");
        }
        const auto &helix_scan_traces =
            m_cfg.whiteboard()->template get<std::vector<intersection_trace_t>>(
                "helix_scan");

        std::size_t n_test_tracks{
            std::min(m_cfg.n_tracks(), helix_scan_traces.size())};
        std::cout << "INFO: Running helix navigation check on: "
                  << m_names.at(0) << "...\n"
                  << std::endl;

        /// Error statistic
        std::size_t n_close_miss{0u}, n_fatal{0u};

        std::ios_base::openmode io_mode = std::ios::trunc | std::ios::out;
        detray::io::file_handle debug_file{"./helix_navigation.txt", io_mode};

        for (auto intersection_trace : helix_scan_traces) {

            if (n_tracks >= m_cfg.n_tracks()) {
                break;
            }

            // Retrieve the test helix
            const auto &start = intersection_trace.front();

            // Follow the test ray with the same track and check, if we find
            // the same volumes and distances along the way
            free_track_parameters_t track(start.track_param);
            detail::helix<algebra_t> helix(track, &B);

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
                // The navigator does not record the initial track position,
                // add it as a dummy record
                obj_tracer.object_trace.insert(
                    obj_tracer.object_trace.begin(),
                    intersection_trace.front().intersection);
                success &=
                    detail::compare_traces(intersection_trace, obj_tracer,
                                           helix, n_tracks, n_test_tracks);
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

            EXPECT_TRUE(success) << "\nFailed on helix " << n_tracks << "/"
                                 << n_test_tracks << ": " << helix << "\n\n";

            ++n_tracks;
        }

        if (n_close_miss > 0u || n_fatal > 0u) {
            std::cout << "-----------------------------------"
                      << "Error Statistic:\n\n"
                      << "\n total: " << n_tracks
                      << "\n (close misses: " << n_close_miss
                      << ", fatal failures: " << n_fatal << ")\n"
                      << "-----------------------------------\n"
                      << std::endl;
        } else {
            std::cout << "-----------------------------------\n"
                      << "Tested " << n_tracks << " helices: OK\n"
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
    /// Truth intersections to compare against
    std::vector<intersection_trace_t> m_intersection_traces;
};

}  // namespace detray::test
