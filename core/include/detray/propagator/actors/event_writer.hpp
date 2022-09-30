/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

struct event_writer : actor {

    struct state {
        state(int event_id, bool write_hits = true,
              bool write_measurements = true)
            : m_write_hits(write_hits),
              m_write_measurements(write_measurements) {
            std::string hit_filename = std::to_string(event_id) + "-hits.csv";
            std::string meas_filename =
                std::to_string(event_id) + "-measurements.csv";

            m_hit_writer = hit_writer(hit_filename);
            m_meas_writer = measurement_writer(meas_filename);
        }

        int particle_id = 0;
        bool m_write_hits = true;
        bool m_write_measurements = true;

        hit_writer m_hit_writer;
        measaurement_writer m_meas_writer;
    };

    struct measurement_kernel {
        using output_type = point2;

        template <typename mask_group_t, typename surface_t,
                  typename bound_params_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const mask_group_t& mask_group, const surface_t& surface,
            bound_params_t& bound_params) {

            const auto& mask = mask_group[surface.mask_range()];

            auto local_coordinate = mask.local();

            return local_coordinate.get_measurement(bound_params);
        }
    }

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& writer_state,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // @note: This should be changed for active surface rather than module
        if (navigation.is_on_module()) {
            if (m_write_hits) {
                csv_hit hit;

                const auto track = stepping();
                const auto pos = track.pos();
                const auto dir = track.dir();
                const auto p = track.p();
                const auto mom = p * dir;
                const auto time = track.time();
                const auto charge = track.charge();

                hit.particle_id = writer_state.particle_id;
                hit.geometry_id = navigation.current()->index;
                hit.tx = pos[0];
                hit.ty = pos[1];
                hit.tz = pos[2];
                hit.tt = time;
                hit.tpx = mom[0];
                hit.tpy = mom[1];
                hit.tpz = mom[2];

                m_hit_writer.append(hit);
            }
            if (m_write_measurements) {
                csv_measurement meas;

                const auto bound_params = stepping._bound_params;
                const auto& det = navigation.detector();
                const auto& mask_store = det->mask_store();
                const auto& is = navigation.current();
                const auto& surface = det->surface_by_index(is->index);

                const auto local =
                    mask_store.template execute<measurement_kernel>(
                        surface.mask_type(), surface, bound_params);

                meas.geometry_id = navigation.current()->index;
                meas.local0 = local[0];
                meas.local1 = local[1];
                meas.phi = bound_params.phi();
                meas.theta = bound_params.theta();
                meas.time = bound_params.time();

                m_meas_writer.append(meas);
            }
        }
    };

}  // namespace detray