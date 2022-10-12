/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/io/utils.hpp"
#include "detray/utils/math.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/simulator.hpp"
#include "tests/common/tools/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

#include <iostream>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

struct global_to_local {

    using output_type = point2;

    template <typename mask_group_t, typename index_t, typename transform_store_t,
              typename surface_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& /*index*/, const transform_store_t& trf_store,
        const surface_t& surface, const point3 pos, const vector3 dir) {

        const auto& trf3 = trf_store[surface.transform()];

        const auto& mask = mask_group[surface.mask_range()];

        auto local_coordinate = mask.local_frame();

        return local_coordinate.global_to_local(trf3, pos, dir);
    }
};

TEST(check_simulation, toy_geometry) {

    // Create geometry
    vecmem::host_memory_resource host_mr;
    const auto detector = create_toy_geometry(host_mr);

    // Create B field
    const vector3 b{0, 0, 2 * unit_constants::T};
    constant_magnetic_field<> b_field(b);

    // Create track generator
    constexpr unsigned int theta_steps{50};
    constexpr unsigned int phi_steps{50};
    const vector3 ori{0, 0, 0};
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, 1 * unit_constants::GeV);

    // Create smearer
    measurement_smearer<scalar> smearer(67 * unit_constants::um,
                                        170 * unit_constants::um);

    std::size_t n_events = 1;
    auto sim = simulator(n_events, detector, b_field, generator, smearer);

    // Do the simulation
    sim.run();

    std::vector<csv_particle> particles;
    std::vector<csv_hit> hits;
    std::vector<csv_measurement> measurements;

    // Check particle data
    const auto io_particles_file = get_event_filename(0, "-particles.csv");
    particle_reader preader(io_particles_file);

    csv_particle io_particle;
    while (preader.read(io_particle)) {
        particles.push_back(io_particle);
    }

    ASSERT_EQ(particles.size(), theta_steps * phi_steps);

    // Check hit & measurement data
    const auto io_hits_file = get_event_filename(0, "-hits.csv");
    hit_reader hreader(
        io_hits_file,
        {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
         "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"});

    csv_hit io_hit;
    while (hreader.read(io_hit)) {
        hits.push_back(io_hit);
    }

    const auto io_measurements_file =
        get_event_filename(0, "-measurements.csv");
    measurement_reader mreader(
        io_measurements_file,
        {"geometry_id", "local_key", "local0", "local1", "phi", "theta", "time",
         "var_local0", "var_local1", "var_phi", "var_theta", "var_time"});

    csv_measurement io_measurement;
    while (mreader.read(io_measurement)) {
        measurements.push_back(io_measurement);
    }

    ASSERT_TRUE(measurements.size() > 0);
    ASSERT_EQ(hits.size(), measurements.size());

    // Let's check if measurement smearing works correctly...
    const auto& trf_store = detector.transform_store();
    const auto& mask_store = detector.mask_store();

    std::vector<scalar> local0_diff;
    std::vector<scalar> local1_diff;

    for (std::size_t i = 0; i < hits.size(); i++) {
        const auto& surface = detector.surface_by_index(hits[i].geometry_id);
        const point3 pos{hits[i].tx, hits[i].ty, hits[i].tz};
        const vector3 mom{hits[i].tpx, hits[i].tpy, hits[i].tpz};
        const auto truth_local = mask_store.template call<global_to_local>(
            surface.mask(), trf_store, surface, pos,
            vector::normalize(mom));

        local0_diff.push_back(truth_local[0] - measurements[i].local0);
        local1_diff.push_back(truth_local[1] - measurements[i].local1);
    }

    const auto var0 = get_variance(local0_diff);
    const auto var1 = get_variance(local1_diff);

    EXPECT_NEAR((std::sqrt(var0) - smearer.stddev[0]) / smearer.stddev[0], 0,
                0.1);
    EXPECT_NEAR((std::sqrt(var1) - smearer.stddev[1]) / smearer.stddev[1], 0,
                0.1);
}