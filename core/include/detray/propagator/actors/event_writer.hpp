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
#include "detray/propagator/base_actor.hpp"

// System include(s).
#include <iomanip>
#include <random>
#include <sstream>

namespace detray {

inline std::string get_event_filename(std::size_t event,
                                      const std::string& suffix) {
    std::stringstream stream;
    stream << "event";
    stream << std::setfill('0') << std::setw(9) << event;
    stream << suffix;
    return stream.str();
}

template <typename scalar_t>
struct measurement_smearer {

    measurement_smearer(scalar_t stddev_local0, scalar_t stddev_local1)
        : stddev({stddev_local0, stddev_local1}) {}

    measurement_smearer(measurement_smearer<scalar_t>& smearer)
        : stddev(smearer.stddev) {}

    std::array<scalar_t, 2> stddev;

    std::random_device rd{};
    std::mt19937 generator{rd()};

    template <int ID>
    scalar_t get() {
        static_assert(ID == 0 || ID == 1);

        if constexpr (ID == 0) {
            return std::normal_distribution<scalar_t>(0, stddev[0])(generator);
        } else if (ID == 1) {
            return std::normal_distribution<scalar_t>(0, stddev[1])(generator);
        }
    }
};

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

        template <typename mask_group_t, typename surface_t>
        inline output_type operator()(
            const mask_group_t& mask_group, const surface_t& surface,
            const bound_track_parameters<transform3_t>& bound_params,
            smearer_t smearer) {

            const auto& mask = mask_group[surface.mask_range()];

            auto local_coordinate = mask.local();

            return local_coordinate.get_measurement(bound_params, smearer);
        }
    };

    template <typename propagator_state_t>
    void operator()(state& writer_state,
                    propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // @note: This should be changed for active surface rather than module
        if (navigation.is_on_module()) {

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
            const auto& det = navigation.detector();
            const auto& mask_store = det->mask_store();
            const auto& is = navigation.current();
            const auto& surface = det->surface_by_index(is->index);

            const auto local = mask_store.template execute<measurement_kernel>(
                surface.mask_type(), surface, bound_params,
                writer_state.m_meas_smearer);

            meas.geometry_id = navigation.current()->index;

            // @todo: apply measurement error
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