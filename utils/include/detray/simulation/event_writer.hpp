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
#include "detray/simulation/measurement_smearer.hpp"

namespace detray {

template <typename transform3_t, typename smearer_t>
struct event_writer : actor {

    using scalar_type = typename transform3_t::scalar_type;

    struct state {
        state(std::size_t event_id, smearer_t& smearer,
              const std::string directory)
            : m_particle_writer(directory +
                                get_event_filename(event_id, "-particles.csv")),
              m_hit_writer(directory +
                           get_event_filename(event_id, "-hits.csv")),
              m_meas_writer(directory +
                            get_event_filename(event_id, "-measurements.csv")),
              m_meas_hit_id_writer(
                  directory +
                  get_event_filename(event_id, "-measurement-simhit-map.csv")),
              m_meas_smearer(smearer) {}

        std::size_t particle_id = -1;
        particle_writer m_particle_writer;
        hit_writer m_hit_writer;
        measurement_writer m_meas_writer;
        meas_hit_id_writer m_meas_hit_id_writer;
        std::size_t m_hit_count = 0;
        smearer_t m_meas_smearer;

        void set_seed(const std::size_t sd) { m_meas_smearer.set_seed(sd); }

        void write_particle(const free_track_parameters<transform3_t>& track) {
            particle_id++;

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
        }
    };

    struct measurement_kernel {

        template <typename mask_group_t, typename index_t>
        inline std::array<scalar_type, 2> operator()(
            const mask_group_t& /*mask_group*/, const index_t& /*index*/,
            const bound_track_parameters<transform3_t>& bound_params,
            smearer_t& smearer) const {

            return smearer(mask_group_t::value_type::shape::name,
                           mask_group_t::value_type::shape::meas_dim,
                           smearer.get_offset(), bound_params);
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

            hit.particle_id = writer_state.particle_id;
            hit.geometry_id = navigation.current_object();
            hit.tx = pos[0];
            hit.ty = pos[1];
            hit.tz = pos[2];
            hit.tt = track.time();
            hit.tpx = mom[0];
            hit.tpy = mom[1];
            hit.tpz = mom[2];

            writer_state.m_hit_writer.append(hit);

            // Write measurements
            csv_measurement meas;

            const auto bound_params = stepping._bound_params;
            auto det = navigation.detector();
            const auto& mask_store = det->mask_store();
            const auto& surface = det->surfaces(hit.geometry_id);

            const auto local = mask_store.template visit<measurement_kernel>(
                surface.mask(), bound_params, writer_state.m_meas_smearer);

            meas.measurement_id = writer_state.m_hit_count;
            meas.geometry_id = hit.geometry_id;
            meas.local_key = "unknown";
            meas.local0 = local[0];
            meas.local1 = local[1];
            auto stddev_0 = writer_state.m_meas_smearer.stddev[0];
            auto stddev_1 = writer_state.m_meas_smearer.stddev[1];
            meas.var_local0 = stddev_0 * stddev_0;
            meas.var_local1 = stddev_1 * stddev_1;
            meas.phi = bound_params.phi();
            meas.theta = bound_params.theta();
            meas.time = bound_params.time();

            writer_state.m_meas_writer.append(meas);

            // Write hit measurement map
            csv_meas_hit_id meas_hit_id;
            meas_hit_id.hit_id = writer_state.m_hit_count;
            meas_hit_id.measurement_id = writer_state.m_hit_count;
            writer_state.m_meas_hit_id_writer.append(meas_hit_id);
            writer_state.m_hit_count++;
        }
    }
};
}  // namespace detray