/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
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
#include "tests/common/tools/event_writer.hpp"
#include "tests/common/tools/random_scatterer.hpp"

namespace detray {

template <typename detector_t, typename track_generator_t, typename smearer_t>
struct simulator {

    using scalar_type = typename detector_t::scalar_type;

    struct config {
        scalar_type pathlimit = std::numeric_limits<scalar>::max();
        scalar_type overstep_tolerance = -10 * detray::unit_constants::um;
        scalar_type step_constraint = 10. * detray::unit_constants::mm;
    };

    using transform3 = typename detector_t::transform3;
    using interactor_t = pointwise_material_interactor<transform3>;
    using bfield_type = typename detector_t::bfield_type;

    using actor_chain_type = actor_chain<
        std::tuple, pathlimit_aborter, parameter_transporter<transform3>,
        interactor_t, random_scatterer<interactor_t>,
        parameter_resetter<transform3>, event_writer<transform3, smearer_t>>;

    using navigator_type = navigator<detector_t>;
    using stepper_type = rk_stepper<typename bfield_type::view_t, transform3,
                                    constrained_step<>>;
    using propagator_type =
        propagator<stepper_type, navigator_type, actor_chain_type>;

    simulator(std::size_t events, const detector_t& det,
              track_generator_t& track_gen, smearer_t& smearer,
              const std::string directory = "")
        : m_events(events),
          m_directory(directory),
          m_detector(std::make_unique<detector_t>(det)),
          m_smearer(smearer) {
        m_track_generator = track_gen;
    }

    config& get_config() { return m_cfg; }

    void run() {
        for (std::size_t event_id = 0; event_id < m_events; event_id++) {
            typename event_writer<transform3, smearer_t>::state writer(
                event_id, m_smearer, m_directory);

            for (auto track : m_track_generator) {

                writer.write_particle(track);

                typename pathlimit_aborter::state aborter{m_cfg.pathlimit};
                typename parameter_transporter<transform3>::state transporter{};
                typename interactor_t::state interactor{};
                typename random_scatterer<interactor_t>::state scatterer(
                    interactor);
                typename parameter_resetter<transform3>::state resetter{};
                typename actor_chain_type::state actor_states =
                    std::tie(aborter, transporter, interactor, scatterer,
                             resetter, writer);

                typename propagator_type::state propagation(
                    track, m_detector->get_bfield(), *m_detector, actor_states);

                propagator_type p({}, {});

                // Set overstep tolerance and stepper constraint
                propagation._stepping().set_overstep_tolerance(
                    m_cfg.overstep_tolerance);
                propagation._stepping.template set_constraint<
                    detray::step::constraint::e_accuracy>(
                    m_cfg.step_constraint);
                p.propagate(propagation);
            }
        }
    }

    private:
    config m_cfg;
    std::size_t m_events = 0;
    std::string m_directory = "";
    std::unique_ptr<detector_t> m_detector;
    track_generator_t m_track_generator;
    smearer_t m_smearer;
};

}  // namespace detray