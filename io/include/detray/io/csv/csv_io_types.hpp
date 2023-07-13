/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"

// DFE include(s).
#include <dfe/dfe_io_dsv.hpp>
#include <dfe/dfe_namedtuple.hpp>

// System include(s).
#include <cstdint>
#include <string>

namespace detray {

struct csv_particle {
    std::uint64_t particle_id = 0u;
    int particle_type = 0;
    int process = 0;
    scalar vx = 0.f;
    scalar vy = 0.f;
    scalar vz = 0.f;
    scalar vt = 0.f;
    scalar px = 0.f;
    scalar py = 0.f;
    scalar pz = 0.f;
    scalar m = 0.f;
    scalar q = 0.f;

    DFE_NAMEDTUPLE(csv_particle, particle_id, particle_type, process, vx, vy,
                   vz, vt, px, py, pz, m, q);
};

using particle_reader = dfe::NamedTupleCsvReader<csv_particle>;
using particle_writer = dfe::NamedTupleCsvWriter<csv_particle>;

struct csv_hit {
    std::uint64_t particle_id = 0u;
    std::uint64_t geometry_id = 0u;
    scalar tx = 0.f;
    scalar ty = 0.f;
    scalar tz = 0.f;
    scalar tt = 0.f;
    scalar tpx = 0.f;
    scalar tpy = 0.f;
    scalar tpz = 0.f;
    scalar te = 0.f;
    scalar deltapx = 0.f;
    scalar deltapy = 0.f;
    scalar deltapz = 0.f;
    scalar deltae = 0.f;
    std::uint64_t index = 0u;

    DFE_NAMEDTUPLE(csv_hit, particle_id, geometry_id, tx, ty, tz, tt, tpx, tpy,
                   tpz, te, deltapx, deltapy, deltapz, deltae, index);
};

using hit_reader = dfe::NamedTupleCsvReader<csv_hit>;
using hit_writer = dfe::NamedTupleCsvWriter<csv_hit>;

struct csv_measurement {

    std::uint64_t measurement_id = 0u;
    std::uint64_t geometry_id = 0u;
    std::string local_key = "";
    scalar local0 = 0.f;
    scalar local1 = 0.f;
    scalar phi = 0.f;
    scalar theta = 0.f;
    scalar time = 0.f;
    scalar var_local0 = 0.f;
    scalar var_local1 = 0.f;
    scalar var_phi = 0.f;
    scalar var_theta = 0.f;
    scalar var_time = 0.f;

    DFE_NAMEDTUPLE(csv_measurement, measurement_id, geometry_id, local_key,
                   local0, local1, phi, theta, time, var_local0, var_local1,
                   var_phi, var_theta, var_time);
};

using measurement_reader = dfe::NamedTupleCsvReader<csv_measurement>;
using measurement_writer = dfe::NamedTupleCsvWriter<csv_measurement>;

struct csv_meas_hit_id {

    std::uint64_t measurement_id = 0u;
    std::uint64_t hit_id = 0u;

    // measurement_id, hit_id
    DFE_NAMEDTUPLE(csv_meas_hit_id, measurement_id, hit_id);
};

using meas_hit_id_reader = dfe::NamedTupleCsvReader<csv_meas_hit_id>;
using meas_hit_id_writer = dfe::NamedTupleCsvWriter<csv_meas_hit_id>;

}  // namespace detray
