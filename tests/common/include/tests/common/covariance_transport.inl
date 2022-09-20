/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/field/constant_magnetic_field.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/resetter.hpp"
#include "detray/propagator/actors/surface_targeter.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <climits>
#include <iostream>

using namespace detray;
using vector2 = __plugin::vector2<scalar>;
using transform3 = __plugin::transform3<scalar>;
using matrix_operator = standard_matrix_operator<scalar>;
using mag_field_t = constant_magnetic_field<>;

constexpr scalar epsilon = 1e-3;

// This tests the covariance transport in rk stepper
TEST(covariance_transport, rk_stepper_cartesian) {

    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    // Use rectangular surfaces
    constexpr bool rectangular = false;

    // Create telescope detector with a single plane
    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {1, 0, 0}, -1};
    std::vector<scalar> positions = {0.,
                                     std::numeric_limits<scalar>::infinity()};

    const auto det = create_telescope_detector<rectangular>(
        host_mr, positions, ln_stepper_t(),
        typename ln_stepper_t::state{default_trk});

    using navigator_t = navigator<decltype(det)>;
    using crk_stepper_t =
        rk_stepper<mag_field_t, transform3, constrained_step<>>;
    using actor_chain_t =
        actor_chain<dtuple, surface_targeter, parameter_transporter<transform3>,
                    resetter<transform3>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain_t>;

    // Generate track starting point
    vector2 local{2, 3};
    vector3 mom{0.001, 0., 0.};
    vector3 dir = vector::normalize(mom);

    scalar time = 0.;
    scalar q = -1.;

    // bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = local[0];
    getter::element(bound_vector, e_bound_loc1, 0) = local[1];
    getter::element(bound_vector, e_bound_phi, 0) = getter::phi(mom);
    getter::element(bound_vector, e_bound_theta, 0) = getter::theta(mom);
    getter::element(bound_vector, e_bound_qoverp, 0) = q / getter::norm(mom);
    getter::element(bound_vector, e_bound_time, 0) = time;

    // bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template zero<e_bound_size, e_bound_size>();
    getter::element(bound_cov, e_bound_loc0, e_bound_loc0) = 1.;
    getter::element(bound_cov, e_bound_loc1, e_bound_loc1) = 1.;
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 1.;
    // Note: Set theta error as ZERO, to constrain the loc1 divergence
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;
    getter::element(bound_cov, e_bound_qoverp, e_bound_qoverp) = 1.;
    getter::element(bound_cov, e_bound_time, e_bound_time) = 1.;

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // B field
    vector3 B{0, 0, 1. * unit_constants::T};
    mag_field_t mag_field(B);

    // Path length per turn
    scalar S = 2. * getter::perp(mom) / getter::norm(B) * M_PI;

    // Actors
    surface_targeter::state targeter{S, 0};
    parameter_transporter<transform3>::state bound_updater{};
    resetter<transform3>::state rst{};

    actor_chain_t::state actor_states = std::tie(targeter, bound_updater, rst);

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, mag_field, det, actor_states);

    crk_stepper_t::state& crk_state = propagation._stepping;

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = 1e-8;

    // RK stepper and its state
    EXPECT_FLOAT_EQ(crk_state().pos()[0], 0);
    EXPECT_FLOAT_EQ(crk_state().pos()[1], local[0]);  // y
    EXPECT_FLOAT_EQ(crk_state().pos()[2], local[1]);  // z
    EXPECT_NEAR(crk_state().dir()[0], dir[0], epsilon);
    EXPECT_NEAR(crk_state().dir()[1], dir[1], epsilon);
    EXPECT_NEAR(crk_state().dir()[2], dir[2], epsilon);

    // helix trajectory
    detail::helix helix(crk_state(), &B);

    // Run propagator
    p.propagate(propagation);

    // RK stepper and its state
    EXPECT_FLOAT_EQ(crk_state().pos()[0], 0);
    EXPECT_NEAR(crk_state().pos()[1], local[0], epsilon);  // y
    EXPECT_NEAR(crk_state().pos()[2], local[1], epsilon);  // z
    EXPECT_NEAR(crk_state().dir()[0], dir[0], epsilon);
    EXPECT_NEAR(crk_state().dir()[1], dir[1], epsilon);
    EXPECT_NEAR(crk_state().dir()[2], dir[2], epsilon);
    EXPECT_NEAR(crk_state.path_length(), targeter._path, epsilon);

    // Bound state after one turn propagation
    const auto bound_param1 = crk_state._bound_params;

    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();

    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (std::size_t i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0),
                    matrix_operator().element(bound_vec1, i, 0), epsilon);
    }

    // covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j),
                        crk_state.path_length() * epsilon);
        }
    }
}

TEST(covariance_transport, linear_stepper_cartesian) {

    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    // Use rectangular surfaces
    constexpr bool unbounded = true;

    // Create telescope detector with a single plane
    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {1, 0, 0}, -1};
    std::vector<scalar> positions = {0., 10., 20., 30., 40., 50., 60.};

    const auto det = create_telescope_detector<unbounded>(
        host_mr, positions, ln_stepper_t(),
        typename ln_stepper_t::state{default_trk});

    using navigator_t = navigator<decltype(det)>;
    using cline_stepper_t = line_stepper<transform3, constrained_step<>>;
    using actor_chain_t = actor_chain<dtuple, parameter_transporter<transform3>,
                                      resetter<transform3>>;
    using propagator_t =
        propagator<cline_stepper_t, navigator_t, actor_chain_t>;

    // Bound vector
    typename bound_track_parameters<transform3>::vector_type bound_vector;
    getter::element(bound_vector, e_bound_loc0, 0) = 0.;
    getter::element(bound_vector, e_bound_loc1, 0) = 0.;
    getter::element(bound_vector, e_bound_phi, 0) = 0.;
    getter::element(bound_vector, e_bound_theta, 0) = M_PI / 4.;
    getter::element(bound_vector, e_bound_qoverp, 0) = -1. / 10.;
    getter::element(bound_vector, e_bound_time, 0) = 0.;

    // Bound covariance
    typename bound_track_parameters<transform3>::covariance_type bound_cov =
        matrix_operator().template identity<e_bound_size, e_bound_size>();

    // Note: Set angle error as ZERO, to constrain the loc0 and loc1 divergence
    getter::element(bound_cov, e_bound_phi, e_bound_phi) = 0.;
    getter::element(bound_cov, e_bound_theta, e_bound_theta) = 0.;

    // Bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    // Actors
    parameter_transporter<transform3>::state bound_updater{};
    resetter<transform3>::state rst{};
    actor_chain_t::state actor_states = std::tie(bound_updater, rst);

    propagator_t p({}, {});
    propagator_t::state propagation(bound_param0, det, actor_states);

    // Run propagator
    p.propagate(propagation);

    // Bound state after one turn propagation
    const auto& bound_param1 = propagation._stepping._bound_params;

    // const auto bound_vec0 = bound_param0.vector();
    // const auto bound_vec1 = bound_param1.vector();

    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the track reaches the final surface
    EXPECT_EQ(bound_param0.surface_link(), 0);
    EXPECT_EQ(bound_param1.surface_link(), 5);

    // Check covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), epsilon);
        }
    }

    // Check covaraince
    for (std::size_t i = 0; i < e_bound_size; i++) {
        for (std::size_t j = 0; j < e_bound_size; j++) {

            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j), epsilon);
        }
    }
}
