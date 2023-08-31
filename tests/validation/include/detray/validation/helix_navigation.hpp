/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/detail/trajectories.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/inspectors.hpp"
#include "detray/validation/detail/navigation_check_helper.hpp"
#include "tests/common/test_base/fixture_base.hpp"
#include "tests/common/tools/particle_gun.hpp"

// System include(s)
#include <iostream>
#include <string>

namespace detray {

/// @brief Test class that runs the helix navigation check on a given detector.
///
/// @note The lifetime of the detector needs to be guaranteed.
template <typename detector_t>
class helix_navigation : public test::fixture_base<> {

    using scalar_t = dscalar<detray::array<typename detector_t::scalar_type>>;
    using transform3_t =
        dtransform3D<detray::array<typename detector_t::scalar_type>>;
    using free_track_parameters_t = free_track_parameters<transform3_t>;

    public:
    using fixture_type = test::fixture_base<>;

    struct config : public fixture_type::configuration {
        using trk_gen_config_t = typename uniform_track_generator<
            free_track_parameters_t>::configuration;

        std::string m_name{"helix_navigation"};
        trk_gen_config_t m_trk_gen_cfg{};

        /// Getters
        /// @{
        const std::string &name() const { return m_name; }
        trk_gen_config_t &track_generator() { return m_trk_gen_cfg; }
        const trk_gen_config_t &track_generator() const {
            return m_trk_gen_cfg;
        }
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
        m_cfg.overstepping_tolerance(cfg.overstepping_tolerance());
        m_cfg.name(cfg.name());
        m_cfg.track_generator() = cfg.track_generator();
    }

    /// Run the check
    void TestBody() override {
        using namespace detray;
        using namespace navigation;

        /// Type that holds the intersection information
        using intersection_t =
            intersection2D<typename detector_t::surface_type, scalar_t, array>;

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
        using b_field_t = typename detector_t::bfield_type;
        using stepper_t = rk_stepper<typename b_field_t::view_t, transform3_t>;
        // Propagator with pathlimit aborter
        using actor_chain_t = actor_chain<dtuple, pathlimit_aborter>;
        using propagator_t = propagator<stepper_t, navigator_t, actor_chain_t>;

        static_assert(
            std::is_same_v<mask<rectangle2D<>>::shape::
                               template intersector_type<intersection_t>,
                           plane_intersector<intersection_t>>,
            "Oops");

        // Propagator
        propagator_t prop(stepper_t{}, navigator_t{});

        // B-field vector for helix
        const auto &origin = m_cfg.track_generator().origin();
        typename b_field_t::view_t bfield_view{m_det.get_bfield()};
        const typename b_field_t::output_t bvec =
            bfield_view.at(origin[0], origin[1], origin[2]);
        const typename fixture_type::point3 B{bvec[0], bvec[1], bvec[2]};

        // Iterate through uniformly distributed momentum directions
        std::size_t n_tracks{0u};
        auto trk_state_generator =
            uniform_track_generator<free_track_parameters_t>(
                m_cfg.track_generator());

        std::cout << "INFO: Running helix navigation check on: "
                  << m_names.at(0) << "\n(" << trk_state_generator.size()
                  << " helices) ...\n"
                  << std::endl;

        for (auto track : trk_state_generator) {

            // Prepare for overstepping in the presence of b fields
            track.set_overstep_tolerance(m_cfg.overstepping_tolerance());

            // Get ground truth helix from track
            detail::helix helix(track, &B);

            // Shoot helix through the detector and record all surfaces it
            // encounters
            const auto intersection_trace = particle_gun::shoot_particle(
                m_det, helix, 15.f * unit<scalar_t>::um);

            // Build actor and propagator states
            pathlimit_aborter::state pathlimit_aborter_state{3.5f *
                                                             unit<scalar_t>::m};
            auto actor_states = std::tie(pathlimit_aborter_state);

            typename propagator_t::state propagation(track, m_det.get_bfield(),
                                                     m_det);

            // Access to navigation information
            auto &inspector = propagation._navigation.inspector();
            auto &obj_tracer = inspector.template get<object_tracer_t>();
            auto &debug_printer = inspector.template get<print_inspector>();

            // Run the propagation
            ASSERT_TRUE(prop.propagate(propagation, actor_states))
                << "\nFailed on helix " << n_tracks << "/"
                << trk_state_generator.size() << ": " << helix << "\n\n"
                << debug_printer.to_string();

            ASSERT_TRUE(detail::compare_traces(intersection_trace, obj_tracer,
                                               helix, n_tracks,
                                               trk_state_generator.size()))
                << debug_printer.to_string();

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
