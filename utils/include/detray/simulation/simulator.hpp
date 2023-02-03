/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/simulation/event_writer.hpp"
#include "detray/simulation/random_scatterer.hpp"

// System include(s).
#include <memory>

namespace detray {

template <typename detector_t, typename track_generator_t, typename smearer_t>
struct simulator {

    using scalar_type = typename detector_t::scalar_type;

    struct config {
        scalar_type overstep_tolerance = -10 * detray::unit<scalar_type>::um;
        scalar_type step_constraint = 10. * detray::unit<scalar_type>::mm;
    };

    using transform3 = typename detector_t::transform3;
    using interactor_t = pointwise_material_interactor<transform3>;
    using bfield_type = typename detector_t::bfield_type;

    using material_actor_t =
        composite_actor<dtuple, interactor_t, random_scatterer<interactor_t>>;

    using actor_chain_type =
        actor_chain<dtuple, parameter_transporter<transform3>, material_actor_t,
                    parameter_resetter<transform3>,
                    event_writer<transform3, smearer_t>>;

    using navigator_type = navigator<detector_t>;
    using stepper_type = rk_stepper<typename bfield_type::view_t, transform3,
                                    constrained_step<>>;
    using propagator_type =
        propagator<stepper_type, navigator_type, actor_chain_type>;

    simulator(unsigned int events, const detector_t& det,
              track_generator_t&& track_gen, smearer_t& smearer,
              const std::string directory = "")
        : m_events(events),
          m_directory(directory),
          m_detector(std::make_unique<detector_t>(det)),
          m_track_generator(
              std::make_unique<track_generator_t>(std::move(track_gen))),
          m_smearer(smearer) {}

    config& get_config() { return m_cfg; }

    void run() {

        for (unsigned int event_id = 0; event_id < m_events; event_id++) {
            typename event_writer<transform3, smearer_t>::state writer(
                event_id, m_smearer, m_directory);

            // Set random seed
            m_scatterer.set_seed(event_id);
            writer.set_seed(event_id);

            auto actor_states = std::tie(m_transporter, m_interactor,
                                         m_scatterer, m_resetter, writer);

            for (auto track : *m_track_generator.get()) {

                writer.write_particle(track);

                typename propagator_type::state propagation(
                    track, m_detector->get_bfield(), *m_detector);

                propagator_type p({}, {});

                // Set overstep tolerance and stepper constraint
                propagation._stepping().set_overstep_tolerance(
                    m_cfg.overstep_tolerance);
                propagation._stepping.template set_constraint<
                    detray::step::constraint::e_accuracy>(
                    m_cfg.step_constraint);

                p.propagate(propagation, actor_states);
            }
        }
    }

    private:
    config m_cfg;
    unsigned int m_events = 0;
    std::string m_directory = "";
    std::unique_ptr<detector_t> m_detector;
    std::unique_ptr<track_generator_t> m_track_generator;
    smearer_t m_smearer;

    /// Actor states
    typename parameter_transporter<transform3>::state m_transporter{};
    typename interactor_t::state m_interactor{};
    typename random_scatterer<interactor_t>::state m_scatterer{};
    typename parameter_resetter<transform3>::state m_resetter{};
};

}  // namespace detray