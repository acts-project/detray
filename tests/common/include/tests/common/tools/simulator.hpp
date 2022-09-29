/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"

namespace detray {

template <typename detector_t, typename field_t, typename track_generator_t>
struct simulator {

    using transform3 = typename detector::tranform3;

    using actor_chain =
        actor_chain<std::tuple, parameter_transporter<transform3>,
                    pointwise_material_interactor<transform3>,
                    parameter_resetter<transform3>>;

    simulator(const detector_t& det, const field_t& field,
              track_generator_t& track_gen)
        : m_detector(det), m_field(field), m_track_generator(track_gen) {}

    void run() {
        for (auto track : m_track_generator) {

            parameter_transporter<transform3>::state transporter{};
            pointwise_material_interactor<transform3>::state interactor{};
            parameter_resetter<transform3>::state resetter{};

            constexpr if (m_track_generator::charge == 0) {}
            else {
            }
        }
    }

    private:
    detector_t m_detector;
    field_t m_field;
    track_generator_t m_track_generator;
};

}  // namespace detray