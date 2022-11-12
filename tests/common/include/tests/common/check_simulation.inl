/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/io/utils.hpp"
#include "detray/utils/statistics.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"
#include "tests/common/tools/create_toy_geometry.hpp"
#include "tests/common/tools/simulator.hpp"
#include "tests/common/tools/track_generators.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <climits>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

struct test_param {
    using point2 = __plugin::point2<scalar>;

    test_param(scalar loc_0, scalar loc_1) {
        loc[0] = loc_0;
        loc[1] = loc_1;
    }

    point2 loc;
    point2 local() const { return loc; }
};

// Test measurement smearer
TEST(check_simulation, measurement_smearer) {

    test_param param(1, 2);
    measurement_smearer<scalar> smearer(0., 0.);

    const auto local_line_1 = smearer("line", 1, {-3, 2}, param);
    ASSERT_FLOAT_EQ(local_line_1[0], 0);
    // Null for one dimensional measurement
    ASSERT_FLOAT_EQ(local_line_1[1], 0);

    const auto local_line_2 = smearer("line", 2, {2, -5}, param);
    ASSERT_FLOAT_EQ(local_line_2[0], 3);
    ASSERT_FLOAT_EQ(local_line_2[1], -3);

    const auto local_rectangle_1 = smearer("rectangle2D", 1, {-3, 2}, param);
    ASSERT_FLOAT_EQ(local_rectangle_1[0], -2);
    // Null for one dimensional measurement
    ASSERT_FLOAT_EQ(local_rectangle_1[1], 0);

    const auto local_rectangle_2 = smearer("rectangle2D", 2, {2, -5}, param);
    ASSERT_FLOAT_EQ(local_rectangle_2[0], 3);
    ASSERT_FLOAT_EQ(local_rectangle_2[1], -3);
}

TEST(check_simulation, toy_geometry) {

    // Create geometry
    vecmem::host_memory_resource host_mr;

    // Create B field
    const vector3 B{0, 0, 2 * unit<scalar>::T};

    // Create geometry
    using b_field_t = decltype(create_toy_geometry(host_mr))::bfield_type;
    const auto detector = create_toy_geometry(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}));

    // Create track generator
    constexpr unsigned int theta_steps{50};
    constexpr unsigned int phi_steps{50};
    const vector3 ori{0, 0, 0};
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, 1 * unit<scalar>::GeV);

    // Create smearer
    measurement_smearer<scalar> smearer(67 * unit<scalar>::um,
                                        170 * unit<scalar>::um);

    std::size_t n_events = 10;
    auto sim = simulator(n_events, detector, generator, smearer);

    // Do the simulation
    sim.run();

    for (std::size_t i_event = 0; i_event < n_events; i_event++) {

        std::vector<csv_particle> particles;
        std::vector<csv_hit> hits;
        std::vector<csv_measurement> measurements;
        std::vector<csv_meas_hit_id> meas_hit_ids;

        // Check particle data
        const auto io_particles_file =
            get_event_filename(i_event, "-particles.csv");
        particle_reader preader(io_particles_file);

        csv_particle io_particle;
        while (preader.read(io_particle)) {
            particles.push_back(io_particle);
        }

        ASSERT_EQ(particles.size(), theta_steps * phi_steps);

        // Check hit & measurement data
        const auto io_hits_file = get_event_filename(i_event, "-hits.csv");
        hit_reader hreader(
            io_hits_file,
            {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
             "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"});

        csv_hit io_hit;
        while (hreader.read(io_hit)) {
            hits.push_back(io_hit);
        }

        const auto io_measurements_file =
            get_event_filename(i_event, "-measurements.csv");
        measurement_reader mreader(
            io_measurements_file,
            {"geometry_id", "local_key", "local0", "local1", "phi", "theta",
             "time", "var_local0", "var_local1", "var_phi", "var_theta",
             "var_time"});

        csv_measurement io_measurement;
        while (mreader.read(io_measurement)) {
            measurements.push_back(io_measurement);
        }

        const auto io_meas_hit_id_file =
            get_event_filename(i_event, "-measurement-simhit-map.csv");
        meas_hit_id_reader meas_hit_id_reader(io_meas_hit_id_file,
                                              {"measurement_id", "hit_id"});
        csv_meas_hit_id io_meas_hit_id;
        while (meas_hit_id_reader.read(io_meas_hit_id)) {
            meas_hit_ids.push_back(io_meas_hit_id);
        }

        ASSERT_TRUE(not measurements.empty());
        ASSERT_EQ(hits.size(), measurements.size());
        ASSERT_EQ(hits.size(), meas_hit_ids.size());

        // Let's check if measurement smearing works correctly...
        std::vector<scalar> local0_diff;
        std::vector<scalar> local1_diff;

        const std::size_t nhits = hits.size();
        for (std::size_t i = 0; i < nhits; i++) {
            const point3 pos{hits[i].tx, hits[i].ty, hits[i].tz};
            const vector3 mom{hits[i].tpx, hits[i].tpy, hits[i].tpz};
            const auto truth_local = detector.global_to_local(
                hits[i].geometry_id, pos, vector::normalize(mom));

            local0_diff.push_back(truth_local[0] - measurements[i].local0);
            local1_diff.push_back(truth_local[1] - measurements[i].local1);

            ASSERT_EQ(meas_hit_ids[i].hit_id, i);
            ASSERT_EQ(meas_hit_ids[i].measurement_id, i);
        }

        const auto var0 = statistics::variance(local0_diff);
        const auto var1 = statistics::variance(local1_diff);

        EXPECT_NEAR((std::sqrt(var0) - smearer.stddev[0]) / smearer.stddev[0],
                    0, 0.1);
        EXPECT_NEAR((std::sqrt(var1) - smearer.stddev[1]) / smearer.stddev[1],
                    0, 0.1);
    }
}

// Test parameters: <initial momentum, theta direction>
class TelescopeDetectorSimulation
    : public ::testing::TestWithParam<std::tuple<scalar, scalar>> {};

TEST_P(TelescopeDetectorSimulation, telescope_detector_simulation) {

    // Create geometry
    vecmem::host_memory_resource host_mr;

    // Build from given module positions
    detail::ray<transform3> traj{{0, 0, 0}, 0, {0, 0, 1}, -1};
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    const auto mat = silicon_tml<scalar>();
    const scalar thickness = 0.17 * unit<scalar>::cm;

    // Detector type
    using detector_type =
        detray::detector<detray::detector_registry::telescope_detector,
                         covfie::field>;

    // Create B field
    const vector3 B{0, 0, 2 * unit<scalar>::T};
    using b_field_t = typename detector_type::bfield_type;

    const auto detector = create_telescope_detector(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        positions, traj, std::numeric_limits<scalar>::infinity(),
        std::numeric_limits<scalar>::infinity(), mat, thickness);

    // Momentum
    const scalar mom = std::get<0>(GetParam());

    // Create track generator
    constexpr unsigned int theta_steps{1};
    constexpr unsigned int phi_steps{1};
    const vector3 ori{0, 0, 0};
    const scalar theta = std::get<1>(GetParam());
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, mom, {theta, theta}, {0., 0.});

    // Create smearer
    measurement_smearer<scalar> smearer(50 * unit<scalar>::um,
                                        50 * unit<scalar>::um);

    std::size_t n_events = 1000;

    auto sim = simulator(n_events, detector, generator, smearer);

    // Run simulation
    sim.run();

    for (std::size_t i_event = 0; i_event < n_events; i_event++) {

        std::vector<csv_measurement> measurements;

        const auto io_measurements_file =
            get_event_filename(i_event, "-measurements.csv");
        measurement_reader mreader(
            io_measurements_file,
            {"geometry_id", "local_key", "local0", "local1", "phi", "theta",
             "time", "var_local0", "var_local1", "var_phi", "var_theta",
             "var_time"});

        csv_measurement io_measurement;
        while (mreader.read(io_measurement)) {
            measurements.push_back(io_measurement);
        }

        // Make sure that number of measurements is equal to the number of
        // physical planes
        ASSERT_EQ(measurements.size(), positions.size() - 2);
    }
}

INSTANTIATE_TEST_SUITE_P(
    Simulation1, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1 * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation2, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation3, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation4, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation5, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1 * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation6, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation7, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation8, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100. * unit<scalar>::GeV, 0.01)));

INSTANTIATE_TEST_SUITE_P(
    Simulation9, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1 * unit<scalar>::GeV, M_PI / 8.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation10, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1. * unit<scalar>::GeV, M_PI / 8.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation11, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10. * unit<scalar>::GeV, M_PI / 8.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation12, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100. * unit<scalar>::GeV, M_PI / 8.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation13, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1 * unit<scalar>::GeV, M_PI / 6.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation14, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1. * unit<scalar>::GeV, M_PI / 6.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation15, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10. * unit<scalar>::GeV, M_PI / 6.)));

INSTANTIATE_TEST_SUITE_P(
    Simulation16, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100. * unit<scalar>::GeV, M_PI / 6.)));