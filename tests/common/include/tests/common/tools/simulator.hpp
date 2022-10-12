/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/event_writer.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/actors/random_scatterer.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

namespace detray {

template <typename detector_t, typename field_t, typename track_generator_t,
          typename smearer_t>
struct simulator {

    using transform3 = typename detector_t::transform3;
    using interactor_t = pointwise_material_interactor<transform3>;

    using actor_chain_type =
        actor_chain<std::tuple, parameter_transporter<transform3>, interactor_t,
                    random_scatterer<interactor_t>,
                    event_writer<transform3, smearer_t>,
                    parameter_resetter<transform3>>;

    using navigator_type = navigator<detector_t>;
    using stepper_type = rk_stepper<field_t, transform3>;
    using propagator_type =
        propagator<stepper_type, navigator_type, actor_chain_type>;

    simulator(std::size_t events, const detector_t& det, const field_t& field,
              track_generator_t& track_gen, smearer_t& smearer)
        : m_events(events),
          m_detector(det),
          m_field(field),
          m_track_generator(track_gen),
          m_smearer(smearer) {}

    void run() {
        for (std::size_t event_id = 0; event_id < m_events; event_id++) {
            typename event_writer<transform3, smearer_t>::state writer(
                event_id, m_smearer);

            for (auto track : m_track_generator) {
                writer.write_particle(track);

                typename parameter_transporter<transform3>::state transporter{};
                typename interactor_t::state interactor{};
                typename random_scatterer<interactor_t>::state scatterer(
                    interactor);
                typename parameter_resetter<transform3>::state resetter{};
                typename actor_chain_type::state actor_states = std::tie(
                    transporter, interactor, scatterer, writer, resetter);

                typename propagator_type::state state(track, m_field,
                                                      m_detector, actor_states);

                propagator_type p({}, {});

                p.propagate(state);
            }
        }
    }

    private:
    std::size_t m_events = 0;
    detector_t m_detector;
    field_t m_field;
    track_generator_t m_track_generator;
    smearer_t m_smearer;
};

}  // namespace detray