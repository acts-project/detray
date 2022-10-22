/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/io/csv_io_types.hpp"
#include "detray/io/utils.hpp"
#include "detray/propagator/base_actor.hpp"
#include "tests/common/tools/measurement_smearer.hpp"

namespace detray {

template <typename transform3_t, typename smearer_t>
struct event_writer : actor {

    struct state {
        state(std::size_t event_id, smearer_t smearer)
            : m_particle_writer(get_event_filename(event_id, "-particles.csv")),
              m_hit_writer(get_event_filename(event_id, "-hits.csv")),
              m_meas_writer(get_event_filename(event_id, "-measurements.csv")),
              m_meas_smearer(smearer) {}

        std::size_t particle_id = 0;
        particle_writer m_particle_writer;
        hit_writer m_hit_writer;
        measurement_writer m_meas_writer;
        smearer_t m_meas_smearer;

        void write_particle(const free_track_parameters<transform3_t>& track) {
            csv_particle particle;
            const auto pos = track.pos();
            const auto mom = track.mom();

            particle.particle_id = particle_id;
            particle.vx = pos[0];
            particle.vy = pos[1];
            particle.vz = pos[2];
            particle.vt = track.time();
            particle.px = mom[0];
            particle.py = mom[1];
            particle.pz = mom[2];
            particle.q = track.charge();

            m_particle_writer.append(particle);
            particle_id++;
        }
    };

    struct measurement_kernel {
        using output_type = point2;

        template <typename mask_group_t, typename index_t, typename surface_t>
        inline output_type operator()(
            const mask_group_t& mask_group, const index_t& /*index*/,
            const surface_t& surface,
            const bound_track_parameters<transform3_t>& bound_params,
            smearer_t smearer) const {

            const auto& mask = mask_group[surface.mask().index()];

            const auto offset = smearer();

            return mask.get_shape().to_measurement(bound_params,
                                                   {offset[0], offset[1]});
        }
    };

    template <typename propagator_state_t>
    void operator()(state& writer_state,
                    propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // triggered only for sensitive surfaces
        if (navigation.is_on_module() &&
            navigation.current()->sf_id == surface_id::e_sensitive) {

            // Write hits
            csv_hit hit;

            const auto track = stepping();
            const auto pos = track.pos();
            const auto mom = track.mom();
            const auto time = track.time();

            hit.particle_id = writer_state.particle_id;
            hit.geometry_id = navigation.current()->index;
            hit.tx = pos[0];
            hit.ty = pos[1];
            hit.tz = pos[2];
            hit.tt = time;
            hit.tpx = mom[0];
            hit.tpy = mom[1];
            hit.tpz = mom[2];

            writer_state.m_hit_writer.append(hit);

            // Write measurements
            csv_measurement meas;

            const auto bound_params = stepping._bound_params;
            auto det = navigation.detector();
            const auto& mask_store = det->mask_store();
            const auto& is = navigation.current();
            const auto& surface = det->surface_by_index(is->index);

            const auto local = mask_store.template call<measurement_kernel>(
                surface.mask(), surface, bound_params,
                writer_state.m_meas_smearer);

            meas.geometry_id = navigation.current()->index;

            meas.local0 = local[0];
            meas.local1 = local[1];
            meas.phi = bound_params.phi();
            meas.theta = bound_params.theta();
            meas.time = bound_params.time();

            writer_state.m_meas_writer.append(meas);
        }
    }
};
}  // namespace detray