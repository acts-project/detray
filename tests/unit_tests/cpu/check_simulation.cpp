/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/detectors/create_telescope_detector.hpp"
#include "detray/detectors/create_toy_geometry.hpp"
#include "detray/io/common/detail/utils.hpp"
#include "detray/masks/masks.hpp"
#include "detray/masks/unbounded.hpp"
#include "detray/simulation/event_generator/track_generators.hpp"
#include "detray/simulation/simulator.hpp"
#include "detray/test/types.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/utils/statistics.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <limits>

using namespace detray;
using transform3 = test::transform3;

constexpr scalar tol{1e-7f};

// Test measurement smearer
GTEST_TEST(detray_simulation, measurement_smearer) {

    const mask<line<false, line_intersector, 1u, true>> ln_1D{
        0u, 10.f * unit<scalar>::mm, 50.f * unit<scalar>::mm};

    const mask<line<false, line_intersector, 2u, true>> ln_2D{
        0u, 10.f * unit<scalar>::mm, 50.f * unit<scalar>::mm};

    const mask<rectangle2D<plane_intersector, 1u, true>> re_1D{
        0u, 10.f * unit<scalar>::mm, 10.f * unit<scalar>::mm};

    const mask<rectangle2D<plane_intersector, 2u, true>> re_2D{
        0u, 10.f * unit<scalar>::mm, 10.f * unit<scalar>::mm};

    bound_track_parameters<transform3> bound_params;
    auto& bound_vec = bound_params.vector();
    getter::element(bound_vec, e_bound_loc0, 0u) = 1.f;
    getter::element(bound_vec, e_bound_loc1, 0u) = 2.f;

    measurement_smearer<transform3> smearer(0.f, 0.f);

    const auto local_line_1 = smearer(ln_1D, {-3.f, 2.f}, bound_params);
    ASSERT_NEAR(local_line_1[0], 0.f, tol);
    // Null for one dimensional measurement
    ASSERT_NEAR(local_line_1[1], 0.f, tol);

    const auto local_line_2 = smearer(ln_2D, {2.f, -5.f}, bound_params);
    ASSERT_NEAR(local_line_2[0], 3.f, tol);
    ASSERT_NEAR(local_line_2[1], -3.f, tol);

    const auto local_rectangle_1 = smearer(re_1D, {-3.f, 2.f}, bound_params);
    ASSERT_NEAR(local_rectangle_1[0], -2.f, tol);
    // Null for one dimensional measurement
    ASSERT_NEAR(local_rectangle_1[1], 0.f, tol);

    const auto local_rectangle_2 = smearer(re_2D, {2.f, -5.f}, bound_params);
    ASSERT_NEAR(local_rectangle_2[0], 3.f, tol);
    ASSERT_NEAR(local_rectangle_2[1], -3.f, tol);
}

GTEST_TEST(detray_simulation, toy_geometry_simulation) {

    // Use deterministic random number generator for testing
    using normal_gen_t =
        random_numbers<scalar, std::normal_distribution<scalar>, std::seed_seq>;

    // Create geometry
    vecmem::host_memory_resource host_mr;
    const auto detector = create_toy_geometry(host_mr);

    // Create track generator
    constexpr unsigned int n_tracks{2500u};
    const vector3 ori{0.f, 0.f, 0.f};
    auto generator =
        random_track_generator<free_track_parameters<transform3>, normal_gen_t>(
            n_tracks, ori);

    // Create smearer
    measurement_smearer<transform3> smearer(67.f * unit<scalar>::um,
                                            170.f * unit<scalar>::um);

    std::size_t n_events{10u};
    auto sim = simulator(n_events, detector, std::move(generator), smearer,
                         test::filenames);

    // Lift step size constraints
    sim.get_config().step_constraint = std::numeric_limits<scalar>::max();

    // Do the simulation
    sim.run();

    for (std::size_t i_event = 0u; i_event < n_events; i_event++) {

        std::vector<csv_particle> particles;
        std::vector<csv_hit> hits;
        std::vector<csv_measurement> measurements;
        std::vector<csv_meas_hit_id> meas_hit_ids;

        // Check particle data
        const auto io_particles_file =
            test::filenames +
            detail::get_event_filename(i_event, "-particles.csv");
        particle_reader preader(io_particles_file);

        csv_particle io_particle;
        while (preader.read(io_particle)) {
            particles.push_back(io_particle);
        }

        ASSERT_EQ(particles.size(), n_tracks);

        // Check hit & measurement data
        const auto io_hits_file =
            test::filenames + detail::get_event_filename(i_event, "-hits.csv");
        hit_reader hreader(
            io_hits_file,
            {"particle_id", "geometry_id", "tx", "ty", "tz", "tt", "tpx", "tpy",
             "tpz", "te", "deltapx", "deltapy", "deltapz", "deltae", "index"});

        csv_hit io_hit;
        while (hreader.read(io_hit)) {
            hits.push_back(io_hit);
        }

        const auto io_measurements_file =
            test::filenames +
            detail::get_event_filename(i_event, "-measurements.csv");
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
            test::filenames +
            detail::get_event_filename(i_event, "-measurement-simhit-map.csv");
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
        for (std::size_t i = 0u; i < nhits; i++) {
            const point3 pos{hits[i].tx, hits[i].ty, hits[i].tz};
            const vector3 mom{hits[i].tpx, hits[i].tpy, hits[i].tpz};
            const auto truth_local =
                detector.global_to_local(geometry::barcode(hits[i].geometry_id),
                                         pos, vector::normalize(mom));

            local0_diff.push_back(truth_local[0] - measurements[i].local0);
            local1_diff.push_back(truth_local[1] - measurements[i].local1);

            ASSERT_EQ(meas_hit_ids[i].hit_id, i);
            ASSERT_EQ(meas_hit_ids[i].measurement_id, i);
        }

        const auto var0 = statistics::variance(local0_diff);
        const auto var1 = statistics::variance(local1_diff);

        EXPECT_NEAR((std::sqrt(var0) - smearer.stddev[0]) / smearer.stddev[0],
                    0.f, 0.1f);
        EXPECT_NEAR((std::sqrt(var1) - smearer.stddev[1]) / smearer.stddev[1],
                    0.f, 0.1f);
    }
}

// Test parameters: <initial momentum, theta direction>
class TelescopeDetectorSimulation
    : public ::testing::TestWithParam<std::tuple<scalar, scalar>> {};

TEST_P(TelescopeDetectorSimulation, telescope_detector_simulation) {

    // Create geometry
    vecmem::host_memory_resource host_mr;

    // Use rectangle surfaces
    mask<rectangle2D<>> rectangle{0u, 1000.f * unit<scalar>::mm,
                                  1000.f * unit<scalar>::mm};

    // Build from given module positions
    std::vector<scalar> positions = {0.f,   50.f,  100.f, 150.f, 200.f, 250.f,
                                     300.f, 350.f, 400.f, 450.f, 500.f};

    const auto mat = silicon_tml<scalar>();
    // A thickness larger than 0.1 cm will flip the track direction of low
    // energy (or non-relativistic) particle due to the large scattering
    const scalar thickness = 0.005f * unit<scalar>::cm;

    // Detector type
    using detector_type = detray::detector<
        detray::detector_registry::template telescope_detector<rectangle2D<>>,
        covfie::field>;

    // Create B field
    const vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};
    using b_field_t = typename detector_type::bfield_type;

    const auto detector = create_telescope_detector(
        host_mr,
        b_field_t(covfie::make_parameter_pack(
            b_field_t::backend_t::configuration_t{B[0], B[1], B[2]})),
        rectangle, positions, mat, thickness);

    // Momentum
    const scalar mom = std::get<0>(GetParam());

    // Create track generator
    constexpr unsigned int theta_steps{1u};
    constexpr unsigned int phi_steps{1u};
    const vector3 ori{0.f, 0.f, 0.f};
    const scalar theta = std::get<1>(GetParam());
    auto generator = uniform_track_generator<free_track_parameters<transform3>>(
        theta_steps, phi_steps, ori, mom, {theta, theta}, {0.f, 0.f});

    // Create smearer
    measurement_smearer<transform3> smearer(50.f * unit<scalar>::um,
                                            50.f * unit<scalar>::um);

    std::size_t n_events{1000u};

    auto sim = simulator(n_events, detector, std::move(generator), smearer,
                         test::filenames);

    // Lift step size constraints
    sim.get_config().step_constraint = std::numeric_limits<scalar>::max();

    // Run simulation
    sim.run();

    for (std::size_t i_event{0u}; i_event < n_events; i_event++) {

        std::vector<csv_measurement> measurements;

        const auto io_measurements_file =
            test::filenames +
            detail::get_event_filename(i_event, "-measurements.csv");
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
        ASSERT_EQ(measurements.size(), positions.size());
    }
}

INSTANTIATE_TEST_SUITE_P(
    Simulation1, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation2, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation3, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation4, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation5, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation6, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation7, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation8, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100.f * unit<scalar>::GeV, 0.01f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation9, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 8.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation10, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 8.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation11, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 8.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation12, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 8.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation13, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(0.1f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 6.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation14, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(1.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 6.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation15, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(10.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 6.f)));

INSTANTIATE_TEST_SUITE_P(
    Simulation16, TelescopeDetectorSimulation,
    ::testing::Values(std::make_tuple(100.f * unit<scalar>::GeV,
                                      constant<scalar>::pi / 6.f)));
