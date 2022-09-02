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
#include "detray/propagator/actors/bound_to_bound_updater.hpp"
#include "detray/propagator/actors/resetter.hpp"
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

using namespace detray;
using vector2 = __plugin::vector2<scalar>;
using vector3 = __plugin::vector3<scalar>;
using point3 = __plugin::point3<scalar>;
using transform3 = __plugin::transform3<scalar>;
using matrix_operator = standard_matrix_operator<scalar>;
using mag_field_t = constant_magnetic_field<>;

constexpr scalar epsilon = 1e-3;
constexpr scalar path_limit = 100 * unit_constants::cm;

// This tests the covariance transport in rk stepper
TEST(covariance_transport, cartesian) {

    vecmem::host_memory_resource host_mr;

    using ln_stepper_t = line_stepper<transform3>;

    // Use rectangular surfaces
    constexpr bool rectangular = false;

    // Create telescope detector with a single plane
    typename ln_stepper_t::free_track_parameters_type default_trk{
        {0, 0, 0}, 0, {1, 0, 0}, -1};
    std::vector<scalar> positions = {0.};

    const auto det = create_telescope_detector<rectangular>(
        host_mr, positions, ln_stepper_t(),
        typename ln_stepper_t::state{default_trk});

    using navigator_t = navigator<decltype(det)>;
    using crk_stepper_t =
        rk_stepper<mag_field_t, transform3, constrained_step<>>;
    using propagator_t = propagator<crk_stepper_t, navigator_t, actor_chain<>>;

    // Generate track starting point
    vector2 local{2, 3};
    vector3 mom{0.02, 0., 0.};
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

    // B field
    vector3 B{0, 0, 1. * unit_constants::T};
    mag_field_t mag_field(B);

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    typename propagator_t::state propagation{bound_param0, mag_field, det};
    crk_stepper_t::state &crk_state = propagation._stepping;
    navigator_t::state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = 1e-8;

    // RK stepper and its state
    crk_stepper_t crk_stepper;

    /*
    ASSERT_FLOAT_EQ(crk_state().pos()[0], 0);
    ASSERT_FLOAT_EQ(crk_state().pos()[1], local[0]);
    ASSERT_FLOAT_EQ(crk_state().pos()[2], local[1]);
    ASSERT_NEAR(crk_state().dir()[0], 1, epsilon);
    ASSERT_NEAR(crk_state().dir()[1], 0, epsilon);
    ASSERT_NEAR(crk_state().dir()[2], 0, epsilon);
    */
    /*

    // test surface
    const vector3 u{0, 1, 0};
    const vector3 w{1, 0, 0};
    const vector3 t{0, 0, 0};
    const transform3 trf(t, w, u);

    // Generate track starting point
    vector3 local{2, 3, 0};
    vector3 mom{0.02, 0., 0.};
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

    // B field
    vector3 B{0, 0, 1. * unit_constants::T};
    mag_field_t mag_field(B);

    // bound track parameter
    const bound_track_parameters<transform3> bound_param0(0, bound_vector,
                                                          bound_cov);

    prop_state<crk_stepper_t::state, nav_state> propagation{
        crk_stepper_t::state(bound_param0, trf, mag_field), nav_state{}};
    crk_stepper_t::state &crk_state = propagation._stepping;
    nav_state &n_state = propagation._navigation;

    // Decrease tolerance down to 1e-8
    crk_state._tolerance = 1e-8;

    // RK stepper and its state
    crk_stepper_t crk_stepper;

    ASSERT_FLOAT_EQ(crk_state().pos()[0], 0);
    ASSERT_FLOAT_EQ(crk_state().pos()[1], local[0]);
    ASSERT_FLOAT_EQ(crk_state().pos()[2], local[1]);
    ASSERT_NEAR(crk_state().dir()[0], 1, epsilon);
    ASSERT_NEAR(crk_state().dir()[1], 0, epsilon);
    ASSERT_NEAR(crk_state().dir()[2], 0, epsilon);

    // helix trajectory
    detail::helix helix(crk_state(), &B);

    // Path length per turn
    scalar S = 2. * getter::norm(mom) / getter::norm(B) * M_PI;

    // Run stepper for one turn
    unsigned int max_steps = 1e4;
    for (unsigned int i = 0; i < max_steps; i++) {

        crk_state.set_constraint(S - crk_state.path_length());

        n_state._step_size = S;

        crk_stepper.step(propagation);

        if (std::abs(S - crk_state.path_length()) < 1e-6) {
            break;
        }

        // Make sure that we didn't reach the end of for loop
        ASSERT_TRUE(i < max_steps - 1);
    }

    // Transport jacobian check

    auto jac_transport = crk_state._jac_transport;
    auto true_J = helix.jacobian(crk_state.path_length());

    for (std::size_t i = 0; i < e_free_size; i++) {
        for (std::size_t j = 0; j < e_free_size; j++) {
            EXPECT_NEAR(matrix_operator().element(jac_transport, i, j),
                        matrix_operator().element(true_J, i, j),
                        crk_state.path_length() * epsilon);
        }
    }

    // Bound parameters check

    // Bound state after one turn propagation
    const auto bound_param1 = crk_stepper.bound_state(propagation, trf);

    const auto bound_vec0 = bound_param0.vector();
    const auto bound_vec1 = bound_param1.vector();

    const auto bound_cov0 = bound_param0.covariance();
    const auto bound_cov1 = bound_param1.covariance();

    // Check if the bound state stays the same after one turn propagation

    // vector
    for (size_type i = 0; i < e_bound_size; i++) {
        EXPECT_NEAR(matrix_operator().element(bound_vec0, i, 0),
                    matrix_operator().element(bound_vec1, i, 0), epsilon);
    }

    // covaraince
    for (size_type i = 0; i < e_bound_size; i++) {
        for (size_type j = 0; j < e_bound_size; j++) {
            EXPECT_NEAR(matrix_operator().element(bound_cov0, i, j),
                        matrix_operator().element(bound_cov1, i, j),
                        crk_state.path_length() * epsilon);
        }
    }
    */
}